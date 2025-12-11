import logging
from pathlib import Path

import pandas as pd

from .io_utils import load_marker_tsv, sanitize_allele, write_oa, write_json, load_sample_sheet, load_oa_table
from .qc import (
    QCReporter,
    add_run_args,
    add_marker_info,
    init_rs_rows,
    update_rs_with_apt,
    add_imputed_section,
    add_imputed_hits_detail,
    write_rs_qc,
    compare_with_genocare,
)
from .vcf_parser import VCFParser


def run_extract(conf: dict, args, logger: logging.Logger) -> int:
    product = args.product
    apt_vcf = conf["apt_vcf"]
    imputed_vcf = conf["imputed_vcf"]
    rsmarker_list = conf["rsmarker_list"]
    sample_sheet_path = args.sample_sheet
    output_dir = Path(conf["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    output_oa_file = output_dir / f"{product}.oa.txt"
    output_json_file = output_oa_file.with_suffix(".json")
    qc_report_path = conf.get("qc_report") or f"{output_oa_file}.qc.report.txt"

    qc = QCReporter(logger)
    add_run_args(qc, product, args)

    sample_sheet = load_sample_sheet(sample_sheet_path, logger)
    sample_names = sample_sheet["sample_names"]
    plate_barcode = sample_sheet["plate_barcode"]
    sample_id_map = sample_sheet["sample_id_map"]

    sample_analysis = qc.analyze_list_file(sample_names, "Sample Sheet (SAMPLE_NAME)")
    sample_list = sample_analysis["unique_items"]
    if len(sample_list) != len(sample_names):
        logger.error("Duplicate SAMPLE_NAME entries detected in sample sheet; aborting.")
        return 1

    marker_df = load_marker_tsv(rsmarker_list, logger)
    marker_list = marker_df["RS_ID"].tolist()
    qc.analyze_list_file(marker_list, "RS Marker List")
    qc_rs_rows = init_rs_rows(marker_list)

    add_marker_info(qc, marker_df)

    qc.add_section("APT VCF 파싱")
    qc.add_key_value("파일", apt_vcf)
    apt_parser = VCFParser(apt_vcf, logger)
    apt_parser.parse()
    qc.add_key_value("VCF 내 샘플 수", len(apt_parser.get_samples()))
    qc.add_key_value("VCF 내 variant 수", len(apt_parser.data))

    qc.add_section("APT VCF에서 마커 추출")
    apt_results, marker_positions = apt_parser.extract_markers(marker_df, sample_list)
    qc.add_key_value("요청 마커 수", len(marker_list))
    qc.add_key_value("APT VCF에서 찾은 마커 수", len(apt_results))
    update_rs_with_apt(qc_rs_rows, marker_list, apt_results, marker_positions)
    missing_in_apt = [m for m in marker_list if m not in apt_results]
    if missing_in_apt:
        qc.add_line("")
        qc.add_line(f"  [APT VCF에 없는 마커: {len(missing_in_apt)}개]")
        for m in missing_in_apt:
            qc.add_line(f"    - {m}")
        logger.info(f"APT VCF에 없는 마커 목록({len(missing_in_apt)}개): {', '.join(missing_in_apt)}")

    imputed_parser = VCFParser(imputed_vcf, logger)
    imputed_parser.parse_header_only()
    imputed_samples = imputed_parser.get_samples()

    needed_positions = set()
    for _, samples_data in apt_results.items():
        for _, geno_info in samples_data.items():
            if geno_info["is_no_call"]:
                chrom = str(geno_info["chrom"]).replace("chr", "")
                pos = str(geno_info["pos"])
                needed_positions.add((chrom, pos))

    missing_positions = {}
    if {"CHROM", "POS"}.issubset(set(marker_df.columns)):
        for _, row in marker_df[marker_df["RS_ID"].isin(missing_in_apt)].iterrows():
            chrom = str(row["CHROM"]).replace("chr", "")
            pos = str(row["POS"])
            missing_positions[row["RS_ID"]] = (chrom, pos)

    missing_position_set = set(missing_positions.values())
    total_positions = needed_positions | missing_position_set

    add_imputed_section(qc, imputed_vcf, imputed_parser, len(total_positions), getattr(args, "imputed_mode", "auto"))
    if total_positions:
        qc.add_line("")
        qc.add_line("  [Imputed 조회 위치 목록 (chrom:pos)]")
        for chrom, pos in sorted(
            total_positions, key=lambda x: (str(x[0]), int(x[1]) if str(x[1]).isdigit() else str(x[1]))
        ):
            qc.add_line(f"    - {chrom}:{pos}")
        logger.info(
            "Imputed 조회 위치 목록: "
            + ", ".join([f"{c}:{p}" for c, p in sorted(total_positions, key=lambda x: (str(x[0]), str(x[1])))])
        )

    imputed_parser.index_positions(
        total_positions,
        target_samples=set(sample_list),
        missing_rs_ids=set(),
        mode=getattr(args, "imputed_mode", "auto"),
    )
    add_imputed_hits_detail(qc, imputed_parser)

    qc.add_section("No Call 처리")
    no_call_count = 0
    filled_from_imputed = 0
    filled_from_normal = 0

    for rs_id, samples_data in apt_results.items():
        for sample_id, geno_info in samples_data.items():
            if geno_info["is_no_call"]:
                no_call_count += 1
                chrom = geno_info["chrom"]
                pos = geno_info["pos"]
                imputed_variant = imputed_parser.get_variant_by_position(chrom, pos)
                if imputed_variant is not None and sample_id in imputed_samples:
                    gt = imputed_variant[sample_id]
                    ref = imputed_variant["REF"]
                    alt = imputed_variant["ALT"]
                    a1, a2 = imputed_parser.decode_genotype(gt, ref, alt)
                    if a1 is not None and a2 is not None:
                        geno_info["allele1"] = a1
                        geno_info["allele2"] = a2
                        geno_info["is_no_call"] = False
                        geno_info["source"] = "IMPUTED"
                        filled_from_imputed += 1
                        continue
                marker_info = marker_df[marker_df["RS_ID"] == rs_id]
                if len(marker_info) > 0:
                    normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                    geno_info["allele1"] = normal_allele
                    geno_info["allele2"] = normal_allele
                    geno_info["is_no_call"] = False
                    geno_info["source"] = "NORMAL_FILLED"
                    filled_from_normal += 1

    imputed_filled_missing = 0
    for rs_id in missing_in_apt:
        pos_info = missing_positions.get(rs_id)
        imputed_variant = None
        if pos_info:
            chrom, pos = pos_info
            imputed_variant = imputed_parser.get_variant_by_position(chrom, pos)
        if imputed_variant is not None:
            chrom = imputed_variant["CHROM"]
            pos = imputed_variant["POS"]
            ref = imputed_variant["REF"]
            alt = imputed_variant["ALT"]
            apt_results[rs_id] = {}
            for sample_id in sample_list:
                if sample_id in imputed_samples:
                    gt = imputed_variant.get(sample_id, "./.")
                    a1, a2 = imputed_parser.decode_genotype(gt, ref, alt)
                    if a1 is not None and a2 is not None:
                        apt_results[rs_id][sample_id] = {
                            "allele1": a1,
                            "allele2": a2,
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "gt_raw": gt,
                            "is_no_call": False,
                            "source": "IMPUTED",
                        }
                        imputed_filled_missing += 1
                        continue
                marker_info = marker_df[marker_df["RS_ID"] == rs_id]
                if len(marker_info) > 0:
                    normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                    apt_results[rs_id][sample_id] = {
                        "allele1": normal_allele,
                        "allele2": normal_allele,
                        "chrom": ".",
                        "pos": ".",
                        "ref": normal_allele,
                        "alt": ".",
                        "gt_raw": "FILLED",
                        "is_no_call": False,
                        "source": "NORMAL_FILLED",
                    }
        else:
            marker_info = marker_df[marker_df["RS_ID"] == rs_id]
            if len(marker_info) > 0:
                normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                apt_results[rs_id] = {}
                for sample_id in sample_list:
                    apt_results[rs_id][sample_id] = {
                        "allele1": normal_allele,
                        "allele2": normal_allele,
                        "chrom": ".",
                        "pos": ".",
                        "ref": normal_allele,
                        "alt": ".",
                        "gt_raw": "FILLED",
                        "is_no_call": False,
                        "source": "NORMAL_FILLED",
                    }

    qc.add_key_value("총 No Call 수", no_call_count)
    qc.add_key_value("Imputed VCF에서 채운 수", filled_from_imputed)
    qc.add_key_value("NORMAL_ALLELE로 채운 수", filled_from_normal)
    qc.add_key_value("Imputed로 채운 누락 RS 수", imputed_filled_missing)

    for rs_id in missing_in_apt:
        if rs_id not in apt_results:
            marker_info = marker_df[marker_df["RS_ID"] == rs_id]
            if len(marker_info) > 0:
                normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                apt_results[rs_id] = {}
                for sample_id in sample_list:
                    apt_results[rs_id][sample_id] = {
                        "allele1": normal_allele,
                        "allele2": normal_allele,
                        "chrom": ".",
                        "pos": ".",
                        "ref": normal_allele,
                        "alt": ".",
                        "gt_raw": "FILLED",
                        "is_no_call": False,
                        "source": "NORMAL_FILLED",
                    }

    # RS 단위 QC 결과 생성
    def _resolve_pos(rs_id: str):
        if rs_id in marker_positions:
            chrom, pos, _, _ = marker_positions[rs_id]
            return (str(chrom).replace("chr", ""), str(pos))
        if rs_id in missing_positions:
            chrom, pos = missing_positions[rs_id]
            return (str(chrom).replace("chr", ""), str(pos))
        if rs_id in apt_results and apt_results[rs_id]:
            sample_any = next(iter(apt_results[rs_id].values()))
            return (str(sample_any.get("chrom", ".")).replace("chr", ""), str(sample_any.get("pos", ".")))
        return None

    for rs_id in marker_list:
        rs_entry = qc_rs_rows[rs_id]
        pos_tuple = _resolve_pos(rs_id)
        rs_entry["imputed_tabix_hit"] = pos_tuple in imputed_parser.tabix_hits if pos_tuple else False
        rs_entry["imputed_stream_hit"] = pos_tuple in imputed_parser.stream_hits if pos_tuple else False

        sources = set()
        if rs_id in apt_results:
            for sample_name, geno_info in apt_results[rs_id].items():
                sources.add(geno_info.get("source", "APT"))
        rs_entry["used_imputed"] = "IMPUTED" in sources
        rs_entry["used_normal_fill"] = "NORMAL_FILLED" in sources
        if "IMPUTED" in sources:
            rs_entry["final_source"] = "IMPUTED"
        elif "NORMAL_FILLED" in sources:
            rs_entry["final_source"] = "NORMAL_FILLED"
        elif sources:
            rs_entry["final_source"] = "APT"
        else:
            rs_entry["final_source"] = "MISSING"

    qc.add_section("OA 형식 출력")
    oa_rows = []
    for rs_id in marker_list:
        if rs_id in apt_results:
            for sample_name in sample_list:
                if sample_name in apt_results[rs_id]:
                    geno_info = apt_results[rs_id][sample_name]
                    source = geno_info.get("source", "APT")
                    sample_id_value = sample_id_map.get(sample_name, sample_name)
                    oa_row = {
                        "Sample ID": sample_id_value,
                        "Plate Barcode": plate_barcode,
                        "Gene Symbol": ".", # TODO: 추후 수정
                        "NCBI SNP Reference": rs_id,
                        "Assay Name or ID": source,
                        "Allele 1 Call": geno_info["allele1"],
                        "Allele 2 Call": geno_info["allele2"],
                    }
                    oa_rows.append(oa_row)

    result_df = pd.DataFrame(oa_rows)
    sample_id_order = [sample_id_map[s] for s in sample_list]
    result_df["Sample ID"] = pd.Categorical(result_df["Sample ID"], categories=sample_id_order, ordered=True)
    result_df["NCBI SNP Reference"] = pd.Categorical(
        result_df["NCBI SNP Reference"], categories=marker_list, ordered=True
    )
    result_df = result_df.sort_values(["Sample ID", "NCBI SNP Reference"]).reset_index(drop=True)

    json_rows = []
    for rs_id in marker_list:
        if rs_id not in apt_results:
            continue
        for sample_name in sample_list:
            if sample_name not in apt_results[rs_id]:
                continue
            geno_info = apt_results[rs_id][sample_name]
            a1 = sanitize_allele(geno_info["allele1"])
            a2 = sanitize_allele(geno_info["allele2"])
            json_rows.append(
                {
                    "sampleId": sample_id_map.get(sample_name, sample_name),
                    "plateBarcode": plate_barcode,
                    "itemCd": rs_id,
                    "result1": a1,
                    "result2": a2,
                }
            )

    qc.add_key_value("총 출력 row 수", len(result_df))
    qc.add_key_value("샘플 수", len(sample_list))
    qc.add_key_value("마커 수", len(marker_list))

    output_path = output_oa_file
    write_oa(
        output_path,
        meta={
            "Source": f"{product} Markers Parser",
            "Mode": product,
            "APT VCF": apt_vcf,
            "Imputed VCF": imputed_vcf,
            "Marker List": rsmarker_list,
            "Sample Sheet": sample_sheet_path,
            "Plate Barcode": plate_barcode,
            "Generated Date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        },
        data=result_df,
        logger=logger,
    )
    write_json(output_json_file, json_rows, logger)

    # QC RS 상태 TSV 저장
    qc_rs_df, qc_rs_path = write_rs_qc(output_dir, product, qc_rs_rows, logger)

    # genocare OA와 비교 (옵션)
    if conf.get("genocare_all_gt"):
        genocare_path = Path(conf["genocare_all_gt"])
        if genocare_path.exists():
            compare_df, compare_summary = compare_with_genocare(
                apt_results,
                sample_list,
                marker_list,
                sample_id_map,
                genocare_path,
                logger,
            )
            compare_path = output_dir / f"{product}.genocare_compare.tsv"
            compare_df.to_csv(compare_path, sep="\t", index=False)
            logger.info(f"Written genocare comparison TSV to {compare_path}")

            qc.add_section("Genocare 비교")
            qc.add_key_value("비교 총 건수", compare_summary["total"])
            qc.add_key_value("일치 수", compare_summary["match"])
            qc.add_key_value("불일치 수", compare_summary["mismatch"])
            qc.add_key_value("genocare에 없는 건수", compare_summary["missing_genocare"])
            qc.add_key_value("비교 파일", str(compare_path))
        else:
            logger.warning(f"Genocare OA file not found: {genocare_path}")

    # QC 리포트에 요약 추가
    qc.add_section("Imputed 조회 상세")
    qc.add_key_value("tabix 조회 성공 위치 수", len(imputed_parser.tabix_hits))
    qc.add_key_value("스트리밍 조회 성공 위치 수", len(imputed_parser.stream_hits))
    qc.add_key_value("고유 위치 수", len(set(imputed_parser.tabix_hits) | set(imputed_parser.stream_hits)))
    qc.add_key_value("총 카운트(indexed_variant_count)", imputed_parser.indexed_variant_count)
    qc.add_section("RS 상태 요약")
    qc.add_key_value("IMPUTED 사용 RS 수", int((qc_rs_df["final_source"] == "IMPUTED").sum()))
    qc.add_key_value("NORMAL_FILLED 사용 RS 수", int((qc_rs_df["final_source"] == "NORMAL_FILLED").sum()))
    qc.add_key_value("APT 사용 RS 수", int((qc_rs_df["final_source"] == "APT").sum()))
    qc.add_key_value("MISSING RS 수", int((qc_rs_df["final_source"] == "MISSING").sum()))

    qc.add_section("처리 완료")
    qc.add_key_value("출력 파일", str(output_oa_file))
    qc.add_key_value("JSON 파일", str(output_json_file))
    qc.add_key_value("상태", "SUCCESS")
    qc.print_report()
    qc.write_report(qc_report_path)
    return 0

