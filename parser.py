#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thin wrapper to call atgcu_ichms CLI.
Usage:
    python parser.py extract --product <product> --sample-list <path> --plate-barcode <str>
    python parser.py api-to-lis --product <product> --json <path>
"""
import sys
from atgcu_ichms.cli import main

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Markers Parser

Description:
    This script parses and processes marker data.

Usage:
    python parser.py --mode <mode> --sample-list <sample_list> --plate-barcode <plate_barcode>

Author: Seung-il Yoo
Date: 2025-12-10
"""

import argparse
import logging
import sys
import json
from pathlib import Path
from typing import Tuple, Dict, List, Optional

import pandas as pd


def load_list_file(file_path: str, logger: logging.Logger = None) -> List[str]:
    """
    리스트 파일을 읽어서 리스트로 반환 (한 줄에 하나씩)

    Args:
        file_path: 리스트 파일 경로
        logger: Logger 객체

    Returns:
        문자열 리스트
    """
    logger = logger or logging.getLogger(__name__)
    items = []

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                items.append(line)

    logger.info(f"Loaded {len(items)} items from {file_path}")
    return items


def load_marker_tsv(file_path: str, logger: logging.Logger = None) -> pd.DataFrame:
    """
    TSV 형식의 RS 마커 리스트 파일을 DataFrame으로 로드

    파일 형식:
        RS_ID	NORMAL_ALLELE	RISK_ALLELE
        rs33972313	C	T
        rs28366003	A	G

    Args:
        file_path: TSV 파일 경로
        logger: Logger 객체

    Returns:
        마커 정보 DataFrame (RS_ID, NORMAL_ALLELE, RISK_ALLELE)
    """
    logger = logger or logging.getLogger(__name__)

    df = pd.read_csv(file_path, sep='\t', encoding='utf-8')

    # 필수 컬럼 확인
    required_cols = ['RS_ID', 'NORMAL_ALLELE', 'RISK_ALLELE']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in marker TSV: {missing_cols}")

    logger.info(f"Loaded {len(df)} markers from {file_path}")
    logger.debug(f"Marker columns: {list(df.columns)}")

    return df


def load_conf(conf_path: str, mode: str) -> Dict[str, any]:
    """
    설정 JSON을 로드하고 mode에 해당하는 설정 반환
    """
    conf_path = Path(conf_path)
    if not conf_path.exists():
        raise FileNotFoundError(f"Config file not found: {conf_path}")

    with open(conf_path, 'r', encoding='utf-8') as f:
        conf_obj = json.load(f)

    if mode not in conf_obj:
        raise KeyError(f"Mode '{mode}' not found in config {conf_path}")

    conf = conf_obj[mode]
    required_keys = ["apt_vcf", "imputed_vcf", "rsmarker_list", "output_dir"]
    missing = [k for k in required_keys if k not in conf or not conf[k]]
    if missing:
        raise KeyError(f"Missing required keys in mode '{mode}': {missing}")

    return conf


def sanitize_allele(allele: Optional[str]) -> str:
    """
    genotype allele을 출력용으로 정규화
    - None 또는 ["N", "*", "-"] 등은 "."
    """
    if allele is None:
        return "."
    a = str(allele)
    if a in ["N", "*", "-"]:
        return "."
    return a


class VCFParser:
    """
    VCF 파일 파싱 클래스 (apt-vcf, imputed-vcf 지원)
    """

    NO_CALL = './.'

    def __init__(self, file_path: str, logger: logging.Logger = None):
        """
        VCFParser 초기화

        Args:
            file_path: VCF 파일 경로
            logger: Logger 객체
        """
        self.file_path = Path(file_path)
        self.logger = logger or logging.getLogger(__name__)
        self.meta_lines: List[str] = []
        self.header: List[str] = []
        self.samples: List[str] = []
        self.data: pd.DataFrame = pd.DataFrame()
        self.position_index: Dict[Tuple[str, str], Dict[str, Dict[str, str]]] = {}
        self.rsid_index: Dict[str, Dict[str, Dict[str, str]]] = {}
        self.indexed_variant_count: int = 0

        self.logger.info(f"VCFParser initialized with file: {self.file_path}")

    def parse_header_only(self) -> None:
        """
        메타/헤더만 파싱 (데이터 라인 로드 없음)
        """
        if self.header:
            return

        meta_lines = []
        header_line = None

        import gzip
        open_func = gzip.open if str(self.file_path).endswith('.gz') else open
        mode = 'rt' if str(self.file_path).endswith('.gz') else 'r'

        with open_func(self.file_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    meta_lines.append(line)
                elif line.startswith('#CHROM'):
                    header_line = line
                    break

        self.meta_lines = meta_lines

        if header_line:
            self.header = header_line.lstrip('#').split('\t')
            if 'FORMAT' in self.header:
                format_idx = self.header.index('FORMAT')
                self.samples = self.header[format_idx + 1:]

        self.logger.info(f"Parsed header only: {len(self.samples)} samples detected")

    def parse(self) -> pd.DataFrame:
        """
        VCF 파일을 파싱하여 DataFrame 반환

        Returns:
            VCF 데이터 DataFrame
        """
        self.logger.info(f"Parsing VCF file: {self.file_path}")

        if not self.file_path.exists():
            raise FileNotFoundError(f"VCF file not found: {self.file_path}")

        meta_lines = []
        header_line = None
        data_lines = []

        # gzip 압축 파일 지원
        import gzip
        open_func = gzip.open if str(self.file_path).endswith('.gz') else open
        mode = 'rt' if str(self.file_path).endswith('.gz') else 'r'

        with open_func(self.file_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    meta_lines.append(line)
                elif line.startswith('#CHROM'):
                    header_line = line
                elif line:
                    data_lines.append(line)

        self.meta_lines = meta_lines

        # 헤더 파싱
        if header_line:
            self.header = header_line.lstrip('#').split('\t')
            # FORMAT 이후가 샘플 컬럼
            if 'FORMAT' in self.header:
                format_idx = self.header.index('FORMAT')
                self.samples = self.header[format_idx + 1:]

        self.logger.info(f"Found {len(self.samples)} samples in VCF")
        self.logger.info(f"Parsing {len(data_lines)} variant lines...")

        # 데이터 파싱
        if data_lines:
            data_rows = [line.split('\t') for line in data_lines]
            self.data = pd.DataFrame(data_rows, columns=self.header)

            # CHROM 컬럼에서 'chr' prefix 제거
            if 'CHROM' in self.data.columns:
                self.data['CHROM'] = self.data['CHROM'].str.replace('chr', '', regex=False)

        self.logger.info(f"Parsed {len(self.data)} variants")
        return self.data

    def get_samples(self) -> List[str]:
        """샘플 리스트 반환"""
        return self.samples

    def get_variant_by_rsid(self, rsid: str) -> pd.Series:
        """
        RS ID로 variant 조회

        Args:
            rsid: RS ID (예: rs12345)

        Returns:
            해당 variant의 Series (없으면 None)
        """
        # 인덱스 기반 조회
        if rsid in self.rsid_index:
            entry = self.rsid_index[rsid]
            row = {
                'CHROM': entry['CHROM'],
                'POS': entry['POS'],
                'REF': entry['REF'],
                'ALT': entry['ALT'],
                'ID': rsid,
            }
            row.update(entry['samples'])
            return pd.Series(row)

        # DataFrame 기반 (기존 방식)
        if 'ID' in self.data.columns:
            matches = self.data[self.data['ID'] == rsid]
            if len(matches) > 0:
                return matches.iloc[0]
        return None

    def get_variant_by_position(self, chrom: str, pos: str) -> pd.Series:
        """
        CHROM:POS로 variant 조회

        Args:
            chrom: 염색체 (예: "1", "10")
            pos: 위치

        Returns:
            해당 variant의 Series (없으면 None)
        """
        chrom = str(chrom).replace('chr', '')
        key = (chrom, str(pos))

        # 인덱스 기반 조회 (메모리 절약)
        if key in self.position_index:
            entry = self.position_index[key]
            row = {
                'CHROM': chrom,
                'POS': str(pos),
                'REF': entry['REF'],
                'ALT': entry['ALT'],
            }
            row.update(entry['samples'])
            return pd.Series(row)

        # DataFrame 기반 (기존 방식)
        if 'CHROM' in self.data.columns and 'POS' in self.data.columns:
            matches = self.data[(self.data['CHROM'] == chrom) & (self.data['POS'] == str(pos))]
            if len(matches) > 0:
                return matches.iloc[0]
        return None

    def decode_genotype(self, gt: str, ref: str, alt: str) -> Tuple[str, str]:
        """
        VCF genotype을 실제 allele로 변환

        Args:
            gt: Genotype 문자열 (예: "0/0", "0/1", "1/1", "./.")
            ref: Reference allele
            alt: Alternative allele

        Returns:
            (Allele1, Allele2) 튜플
        """
        if gt in ['./.', '.', './.']:
            return (None, None)  # No call

        # FORMAT 필드가 포함된 경우 genotype만 추출
        gt_core = gt.split(':', 1)[0] if gt else gt

        # | 또는 / 로 분리
        separator = '|' if '|' in gt_core else '/'
        parts = gt_core.split(separator)

        if len(parts) != 2:
            return (None, None)

        alleles = [ref, alt]
        try:
            a1_idx = int(parts[0])
            a2_idx = int(parts[1])

            a1 = alleles[a1_idx] if a1_idx < len(alleles) else None
            a2 = alleles[a2_idx] if a2_idx < len(alleles) else None

            return (a1, a2)
        except (ValueError, IndexError):
            return (None, None)

    def index_positions(
        self,
        positions: set,
        target_samples: Optional[set] = None,
        missing_rs_ids: Optional[set] = None
    ) -> None:
        """
        특정 위치만 인덱싱하여 메모리를 아끼는 함수

        Args:
            positions: {(chrom, pos)} 형태의 세트
            target_samples: 조회할 샘플 서브셋 (None이면 전체)
        """
        normalized_positions = {(str(c).replace('chr', ''), str(p)) for c, p in positions}
        missing_rs_ids = set(missing_rs_ids) if missing_rs_ids else set()
        if not normalized_positions and not missing_rs_ids:
            self.logger.info("No positions or RS IDs requested; skipping imputed VCF indexing.")
            return

        self.parse_header_only()
        self.position_index = {}
        self.rsid_index = {}
        self.indexed_variant_count = 0

        target_samples = set(target_samples) if target_samples else set(self.samples)
        sample_to_idx = {s: i for i, s in enumerate(self.samples)}

        # Tabix 인덱스가 있으면 random access 시도
        tbi_path = Path(str(self.file_path) + ".tbi")
        csi_path = Path(str(self.file_path) + ".csi")
        has_index = tbi_path.exists() or csi_path.exists()

        used_tabix = False
        if has_index and normalized_positions:
            try:
                import pysam

                self.logger.info("Using tabix index for targeted fetching")
                vcf = pysam.VariantFile(str(self.file_path))
                for chrom, pos in normalized_positions:
                    for query_chrom in (chrom, f"chr{chrom}"):
                        try:
                            for rec in vcf.fetch(query_chrom, int(pos) - 1, int(pos)):
                                if str(rec.pos) != pos:
                                    continue
                                ref = rec.ref
                                alt = rec.alts[0] if rec.alts else '.'
                                sample_gts = {}
                                for sample in target_samples:
                                    if sample in rec.samples:
                                        gt_tuple = rec.samples[sample].get('GT')
                                        if gt_tuple:
                                            gt_str = '/'.join('.' if a is None else str(a) for a in gt_tuple)
                                        else:
                                            gt_str = './.'
                                        sample_gts[sample] = gt_str
                                self.position_index[(chrom, pos)] = {
                                    'REF': ref,
                                    'ALT': alt,
                                    'samples': sample_gts
                                }
                                self.indexed_variant_count += 1
                                break
                            else:
                                continue
                            break
                        except ValueError:
                            continue
                used_tabix = True
            except ImportError:
                self.logger.warning("pysam not installed; falling back to streaming read.")
            except Exception as e:
                self.logger.warning(f"Tabix fetch failed ({e}); falling back to streaming read.")

        # 스트리밍 fallback (한 번만 스캔, 필요한 위치만 저장)
        import gzip
        open_func = gzip.open if str(self.file_path).endswith('.gz') else open
        mode = 'rt' if str(self.file_path).endswith('.gz') else 'r'

        # streaming: 실행 조건
        need_stream = True
        if used_tabix and not missing_rs_ids and len(self.position_index) >= len(normalized_positions):
            need_stream = False

        if not need_stream:
            return

        self.logger.info(
            f"Streaming VCF to find positions ({len(normalized_positions)}) "
            f"and RS IDs ({len(missing_rs_ids)})..."
        )
        with open_func(self.file_path, mode, encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                raw_chrom = parts[0]
                chrom = raw_chrom.replace('chr', '')
                pos = parts[1]
                rsid = parts[2] if len(parts) > 2 else ''
                key = (chrom, pos)
                hit_position = key in normalized_positions
                hit_rsid = rsid in missing_rs_ids
                if not hit_position and not hit_rsid:
                    continue

                ref = parts[3]
                alt = parts[4].split(',')[0] if parts[4] else '.'
                format_fields = parts[8].split(':')
                sample_fields = parts[9:]

                sample_gts = {}
                for sample in target_samples:
                    idx = sample_to_idx.get(sample)
                    if idx is None or idx >= len(sample_fields):
                        continue
                    sample_value = sample_fields[idx]
                    gt_value = sample_value.split(':', 1)[0] if sample_value else './.'
                    sample_gts[sample] = gt_value

                if hit_position:
                    self.position_index[key] = {
                        'REF': ref,
                        'ALT': alt,
                        'samples': sample_gts
                    }
                    self.indexed_variant_count += 1

                if hit_rsid:
                    self.rsid_index[rsid] = {
                        'CHROM': chrom,
                        'POS': pos,
                        'REF': ref,
                        'ALT': alt,
                        'samples': sample_gts
                    }
                    self.indexed_variant_count += 1

    def extract_markers(
        self,
        marker_df: pd.DataFrame,
        sample_list: List[str]
    ) -> Dict[str, Dict[str, Dict]]:
        """
        마커 리스트와 샘플 리스트에 해당하는 genotype 추출

        Args:
            marker_df: 마커 정보 DataFrame (RS_ID, NORMAL_ALLELE, RISK_ALLELE)
            sample_list: 추출할 샘플 리스트

        Returns:
            {rs_id: {sample_id: {'allele1': X, 'allele2': Y, 'chrom': Z, 'pos': W}}}
        """
        results = {}
        marker_positions = {}  # {rs_id: (chrom, pos, ref, alt)}

        for _, row in marker_df.iterrows():
            rs_id = row['RS_ID']
            variant = self.get_variant_by_rsid(rs_id)

            if variant is None:
                self.logger.warning(f"RS ID not found in VCF: {rs_id}")
                continue

            chrom = variant['CHROM']
            pos = variant['POS']
            ref = variant['REF']
            alt = variant['ALT']

            marker_positions[rs_id] = (chrom, pos, ref, alt)
            results[rs_id] = {}

            for sample in sample_list:
                if sample in self.samples:
                    gt = variant[sample]
                    allele1, allele2 = self.decode_genotype(gt, ref, alt)

                    results[rs_id][sample] = {
                        'allele1': allele1,
                        'allele2': allele2,
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'gt_raw': gt,
                        'is_no_call': allele1 is None or allele2 is None
                    }
                else:
                    self.logger.warning(f"Sample not found in VCF: {sample}")

        self.logger.info(f"Extracted {len(results)} markers for {len(sample_list)} samples")
        return results, marker_positions


def setup_logger(name: str, log_level: str = "INFO", log_file: Optional[str] = None) -> logging.Logger:
    """
    로거 설정 및 반환

    Args:
        name: 로거 이름
        log_level: 로깅 레벨 (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    Returns:
        설정된 Logger 객체
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, log_level.upper()))

    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console handler
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, log_level.upper()))
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # File handler (옵션)
    if log_file and not any(isinstance(h, logging.FileHandler) and getattr(h, 'baseFilename', '') == str(Path(log_file)) for h in logger.handlers):
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(getattr(logging, log_level.upper()))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


class Parser:
    """
    DTC53 마커 데이터를 파싱하는 클래스
    """

    def __init__(self, input_file: str, output_file: str, logger: logging.Logger = None):
        """
        Parser 초기화

        Args:
            input_file: 입력 파일 경로
            output_file: 출력 파일 경로
            logger: Logger 객체
        """
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)
        self.logger = logger or logging.getLogger(__name__)

        self.logger.info(f"Parser initialized with input: {self.input_file}")

    def validate_input(self) -> bool:
        """
        입력 파일 유효성 검사

        Returns:
            유효하면 True, 아니면 False
        """
        if not self.input_file.exists():
            self.logger.error(f"Input file not found: {self.input_file}")
            return False

        if not self.input_file.is_file():
            self.logger.error(f"Input path is not a file: {self.input_file}")
            return False

        self.logger.debug(f"Input file validated: {self.input_file}")
        return True

    def parse(self) -> dict:
        """
        입력 파일을 파싱

        Returns:
            파싱된 데이터 딕셔너리
        """
        self.logger.info("Starting parsing...")

        data = {}

        # TODO: 실제 파싱 로직 구현
        with open(self.input_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # 파싱 로직 추가
                self.logger.debug(f"Processing line {line_num}: {line[:50]}...")

        self.logger.info("Parsing completed")
        return data

    def write_output(self, data: dict) -> None:
        """
        결과를 출력 파일에 저장

        Args:
            data: 저장할 데이터
        """
        self.logger.info(f"Writing output to: {self.output_file}")

        # 출력 디렉토리 생성
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(self.output_file, 'w', encoding='utf-8') as f:
            # TODO: 실제 출력 로직 구현
            pass

        self.logger.info("Output written successfully")

    def run(self) -> int:
        """
        전체 파싱 프로세스 실행

        Returns:
            성공 시 0, 실패 시 1
        """
        try:
            if not self.validate_input():
                return 1

            data = self.parse()
            self.write_output(data)

            self.logger.info("Process completed successfully")
            return 0

        except Exception as e:
            self.logger.exception(f"Error during processing: {e}")
            return 1


class OAFileParser:
    """
    GenoCareAll_genotype.txt (OA File) 파싱 클래스

    파일 형식:
        - 상단: # 으로 시작하는 메타데이터 라인
        - 하단: Sample ID로 시작하는 탭 구분 테이블
    """

    def __init__(self, file_path: str, logger: logging.Logger = None):
        """
        OAFileParser 초기화

        Args:
            file_path: OA 파일 경로
            logger: Logger 객체
        """
        self.file_path = Path(file_path)
        self.logger = logger or logging.getLogger(__name__)
        self.meta: Dict[str, str] = {}
        self.data: pd.DataFrame = pd.DataFrame()

        self.logger.info(f"OAFileParser initialized with file: {self.file_path}")

    def parse(self) -> Tuple[Dict[str, str], pd.DataFrame]:
        """
        OA 파일을 파싱하여 메타데이터와 데이터프레임 반환

        Returns:
            Tuple[Dict, DataFrame]: (메타데이터 딕셔너리, 데이터 DataFrame)
        """
        self.logger.info(f"Parsing OA file: {self.file_path}")

        if not self.file_path.exists():
            self.logger.error(f"File not found: {self.file_path}")
            raise FileNotFoundError(f"File not found: {self.file_path}")

        meta_lines = []
        data_start_line = 0

        # 파일을 읽어서 메타데이터와 데이터 시작 위치 파악
        with open(self.file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f):
                line = line.strip()

                # 빈 라인 스킵
                if not line:
                    continue

                # 메타데이터 라인 (#으로 시작)
                if line.startswith('#'):
                    meta_lines.append(line)
                else:
                    # 데이터 시작 (Sample ID 헤더)
                    data_start_line = line_num
                    break

        # 메타데이터 파싱
        self.meta = self._parse_meta(meta_lines)
        self.logger.info(f"Parsed {len(self.meta)} metadata entries")

        # 데이터 테이블 파싱 (pandas로 읽기)
        self.data = pd.read_csv(
            self.file_path,
            sep='\t',
            skiprows=data_start_line,
            encoding='utf-8'
        )
        self.logger.info(f"Parsed data table: {len(self.data)} rows, {len(self.data.columns)} columns")
        self.logger.debug(f"Columns: {list(self.data.columns)}")

        return self.meta, self.data

    def _parse_meta(self, meta_lines: list) -> Dict[str, str]:
        """
        메타데이터 라인들을 딕셔너리로 파싱

        Args:
            meta_lines: # 으로 시작하는 라인 리스트

        Returns:
            메타데이터 딕셔너리
        """
        meta = {}

        for line in meta_lines:
            # "# Key : Value" 형식 파싱
            line = line.lstrip('#').strip()

            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                meta[key] = value
            else:
                self.logger.warning(f"Invalid meta line format: {line}")

        return meta

    def get_meta(self) -> Dict[str, str]:
        """메타데이터 딕셔너리 반환"""
        return self.meta

    def get_data(self) -> pd.DataFrame:
        """데이터 DataFrame 반환"""
        return self.data

    def get_samples(self) -> list:
        """고유 샘플 ID 리스트 반환"""
        if 'Sample ID' in self.data.columns:
            return self.data['Sample ID'].unique().tolist()
        return []

    def get_markers(self) -> list:
        """고유 RS 마커 리스트 반환"""
        if 'NCBI SNP Reference' in self.data.columns:
            return self.data['NCBI SNP Reference'].unique().tolist()
        return []

    def filter_data(
        self,
        sample_list: List[str] = None,
        marker_list: List[str] = None
    ) -> pd.DataFrame:
        """
        샘플 리스트와 마커 리스트로 DataFrame 필터링

        Args:
            sample_list: 필터링할 샘플 ID 리스트
            marker_list: 필터링할 RS 마커 리스트

        Returns:
            필터링된 DataFrame
        """
        filtered_df = self.data.copy()

        # 샘플 필터링
        if sample_list is not None and len(sample_list) > 0:
            if 'Sample ID' in filtered_df.columns:
                before_count = len(filtered_df)
                filtered_df = filtered_df[filtered_df['Sample ID'].isin(sample_list)]
                self.logger.info(f"Filtered by samples: {before_count} -> {len(filtered_df)} rows")
            else:
                self.logger.warning("'Sample ID' column not found in data")

        # 마커 필터링
        if marker_list is not None and len(marker_list) > 0:
            if 'NCBI SNP Reference' in filtered_df.columns:
                before_count = len(filtered_df)
                filtered_df = filtered_df[filtered_df['NCBI SNP Reference'].isin(marker_list)]
                self.logger.info(f"Filtered by markers: {before_count} -> {len(filtered_df)} rows")
            else:
                self.logger.warning("'NCBI SNP Reference' column not found in data")

        self.logger.info(f"Final filtered data: {len(filtered_df)} rows")
        return filtered_df

    def filter_and_fill_data(
        self,
        sample_list: List[str],
        marker_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        샘플/마커로 필터링하고, 누락된 마커는 NORMAL_ALLELE로 homozygous genotype 생성

        Args:
            sample_list: 필터링할 샘플 ID 리스트
            marker_df: 마커 정보 DataFrame (RS_ID, NORMAL_ALLELE, RISK_ALLELE 컬럼 필수)

        Returns:
            필터링 및 채워진 DataFrame
        """
        marker_list = marker_df['RS_ID'].tolist()

        # 1. 샘플 필터링
        filtered_df = self.data.copy()
        if sample_list is not None and len(sample_list) > 0:
            if 'Sample ID' in filtered_df.columns:
                before_count = len(filtered_df)
                filtered_df = filtered_df[filtered_df['Sample ID'].isin(sample_list)]
                self.logger.info(f"Filtered by samples: {before_count} -> {len(filtered_df)} rows")

        # 2. 마커 필터링 (존재하는 마커만)
        if 'NCBI SNP Reference' in filtered_df.columns:
            before_count = len(filtered_df)
            filtered_df = filtered_df[filtered_df['NCBI SNP Reference'].isin(marker_list)]
            self.logger.info(f"Filtered by markers: {before_count} -> {len(filtered_df)} rows")

        # 3. 누락된 마커 찾기 및 NORMAL_ALLELE로 채우기
        existing_markers_per_sample = {}
        if 'Sample ID' in filtered_df.columns and 'NCBI SNP Reference' in filtered_df.columns:
            for sample in sample_list:
                sample_data = filtered_df[filtered_df['Sample ID'] == sample]
                existing_markers_per_sample[sample] = sample_data['NCBI SNP Reference'].tolist()

        # 누락된 마커에 대해 homozygous NORMAL_ALLELE 데이터 생성
        fill_rows = []
        for sample in sample_list:
            existing_markers = existing_markers_per_sample.get(sample, [])
            missing_markers = [m for m in marker_list if m not in existing_markers]

            if missing_markers:
                self.logger.info(f"Sample {sample}: filling {len(missing_markers)} missing markers with NORMAL_ALLELE")

                for rs_id in missing_markers:
                    marker_info = marker_df[marker_df['RS_ID'] == rs_id]
                    if len(marker_info) > 0:
                        normal_allele = marker_info.iloc[0]['NORMAL_ALLELE']

                        # OA 형식의 row 생성
                        fill_row = {
                            'Sample ID': sample,
                            'Plate Barcode': 'FILLED',
                            'Gene Symbol': '.',
                            'NCBI SNP Reference': rs_id,
                            'Assay Name or ID': rs_id,
                            'Allele 1 Call': normal_allele,
                            'Allele 2 Call': normal_allele
                        }
                        fill_rows.append(fill_row)

        # 채워진 데이터 추가
        if fill_rows:
            fill_df = pd.DataFrame(fill_rows)
            filtered_df = pd.concat([filtered_df, fill_df], ignore_index=True)
            self.logger.info(f"Added {len(fill_rows)} filled rows for missing markers")

        # 샘플 ID와 RS ID로 정렬
        if 'Sample ID' in filtered_df.columns and 'NCBI SNP Reference' in filtered_df.columns:
            filtered_df = filtered_df.sort_values(['Sample ID', 'NCBI SNP Reference']).reset_index(drop=True)

        self.logger.info(f"Final data after fill: {len(filtered_df)} rows")
        return filtered_df

    def write_oa_format(
        self,
        output_path: str,
        data: pd.DataFrame = None,
        meta: Dict[str, str] = None
    ) -> None:
        """
        OA 형식으로 파일 출력 (메타데이터 + 탭 구분 테이블)

        Args:
            output_path: 출력 파일 경로
            data: 출력할 DataFrame (None이면 self.data 사용)
            meta: 출력할 메타데이터 (None이면 self.meta 사용)
        """
        output_path = Path(output_path)
        data = data if data is not None else self.data
        meta = meta if meta is not None else self.meta

        # 출력 디렉토리 생성
        output_path.parent.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Writing OA format to: {output_path}")

        with open(output_path, 'w', encoding='utf-8') as f:
            # 메타데이터 출력
            for key, value in meta.items():
                f.write(f"# {key} : {value}\n")

            # 메타데이터와 테이블 사이 빈 줄 추가
            f.write("\n")

            # 데이터 테이블 출력 (헤더 포함)
            data.to_csv(f, sep='\t', index=False)

        self.logger.info(f"Written {len(data)} rows to {output_path}")


class QCReporter:
    """
    QC(Quality Control) 분석 리포트 생성 클래스
    """

    def __init__(self, logger: logging.Logger = None):
        """
        QCReporter 초기화

        Args:
            logger: Logger 객체
        """
        self.logger = logger or logging.getLogger(__name__)
        self.report_lines: List[str] = []

    def add_section(self, title: str) -> None:
        """리포트에 섹션 추가"""
        self.report_lines.append("")
        self.report_lines.append("=" * 70)
        self.report_lines.append(f"  {title}")
        self.report_lines.append("=" * 70)

    def add_line(self, text: str) -> None:
        """리포트에 라인 추가"""
        self.report_lines.append(text)

    def add_key_value(self, key: str, value) -> None:
        """리포트에 key-value 쌍 추가"""
        self.report_lines.append(f"  - {key}: {value}")

    def analyze_list_file(
        self,
        items: List[str],
        list_name: str
    ) -> Dict[str, any]:
        """
        리스트 파일 분석 (총 개수, 중복 제거 개수, 중복 항목)

        Args:
            items: 분석할 리스트
            list_name: 리스트 이름 (로깅용)

        Returns:
            분석 결과 딕셔너리
        """
        total_count = len(items)
        unique_items = list(dict.fromkeys(items))
        unique_count = len(unique_items)
        duplicates = [item for item in set(items) if items.count(item) > 1]

        self.add_section(f"{list_name} 분석")
        self.add_key_value("총 항목 수", total_count)
        self.add_key_value("중복 제거 후 항목 수", unique_count)
        self.add_key_value("중복 항목 수", total_count - unique_count)

        if duplicates:
            self.add_line(f"  - 중복된 항목들: {duplicates}")

        return {
            "total": total_count,
            "unique": unique_count,
            "duplicates": duplicates,
            "unique_items": unique_items
        }

    def analyze_oa_metadata(self, meta: Dict[str, str]) -> None:
        """
        OA 파일 메타데이터 분석

        Args:
            meta: 메타데이터 딕셔너리
        """
        self.add_section("Input OA File 메타데이터")

        for key, value in meta.items():
            self.add_key_value(key, value if value else "(empty)")

    def analyze_oa_data(self, data: pd.DataFrame) -> Dict[str, any]:
        """
        OA 파일 데이터 분석

        Args:
            data: OA DataFrame

        Returns:
            분석 결과 딕셔너리
        """
        self.add_section("Input OA File 데이터 분석")

        # 총 row 수
        total_rows = len(data)
        self.add_key_value("총 데이터 row 수", total_rows)

        # 샘플 분석
        if 'Sample ID' in data.columns:
            samples_in_data = data['Sample ID'].unique().tolist()
            self.add_key_value("고유 샘플 수", len(samples_in_data))

            # 샘플별 마커 개수
            sample_marker_counts = data.groupby('Sample ID').size().to_dict()
            self.add_line("")
            self.add_line("  [샘플별 RS 마커 개수]")
            for sample, count in sample_marker_counts.items():
                self.add_line(f"    - {sample}: {count} markers")
        else:
            samples_in_data = []
            sample_marker_counts = {}

        # 마커 분석
        if 'NCBI SNP Reference' in data.columns:
            markers_in_data = data['NCBI SNP Reference'].unique().tolist()
            self.add_key_value("고유 RS 마커 수", len(markers_in_data))
        else:
            markers_in_data = []

        return {
            "total_rows": total_rows,
            "samples": samples_in_data,
            "markers": markers_in_data,
            "sample_marker_counts": sample_marker_counts
        }

    def analyze_filter_results(
        self,
        requested_samples: List[str],
        requested_markers: List[str],
        available_samples: List[str],
        available_markers: List[str],
        filtered_df: pd.DataFrame
    ) -> Dict[str, any]:
        """
        필터링 결과 분석 (누락된 샘플/마커 확인)

        Args:
            requested_samples: 요청된 샘플 리스트
            requested_markers: 요청된 마커 리스트
            available_samples: OA 파일에 존재하는 샘플 리스트
            available_markers: OA 파일에 존재하는 마커 리스트
            filtered_df: 필터링된 DataFrame

        Returns:
            분석 결과 딕셔너리
        """
        self.add_section("필터링 결과 분석")

        # 요청했지만 존재하지 않는 샘플
        missing_samples = [s for s in requested_samples if s not in available_samples]
        found_samples = [s for s in requested_samples if s in available_samples]

        self.add_key_value("요청 샘플 수", len(requested_samples))
        self.add_key_value("매칭된 샘플 수", len(found_samples))
        self.add_key_value("누락된 샘플 수", len(missing_samples))

        if missing_samples:
            self.add_line("")
            self.add_line("  [존재하지 않는 샘플 목록]")
            for sample in missing_samples:
                self.add_line(f"    - {sample}")

        # 요청했지만 존재하지 않는 마커
        missing_markers = [m for m in requested_markers if m not in available_markers]
        found_markers = [m for m in requested_markers if m in available_markers]

        self.add_line("")
        self.add_key_value("요청 마커 수", len(requested_markers))
        self.add_key_value("매칭된 마커 수", len(found_markers))
        self.add_key_value("누락된 마커 수", len(missing_markers))

        if missing_markers:
            self.add_line("")
            self.add_line("  [존재하지 않는 RS 마커 목록]")
            for marker in missing_markers:
                self.add_line(f"    - {marker}")

        # 최종 필터링 결과
        self.add_line("")
        self.add_key_value("최종 필터링된 row 수", len(filtered_df))

        if 'Sample ID' in filtered_df.columns:
            final_samples = filtered_df['Sample ID'].unique().tolist()
            self.add_key_value("최종 샘플 수", len(final_samples))

        if 'NCBI SNP Reference' in filtered_df.columns:
            final_markers = filtered_df['NCBI SNP Reference'].unique().tolist()
            self.add_key_value("최종 마커 수", len(final_markers))

        return {
            "missing_samples": missing_samples,
            "found_samples": found_samples,
            "missing_markers": missing_markers,
            "found_markers": found_markers
        }

    def get_report(self) -> str:
        """전체 리포트 문자열 반환"""
        return "\n".join(self.report_lines)

    def write_report(self, output_path: str) -> None:
        """
        리포트를 파일로 저장

        Args:
            output_path: 출력 파일 경로
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'a', encoding='utf-8') as f:
            f.write(self.get_report())
            f.write("\n")

        self.logger.info(f"QC report written to: {output_path}")

    def print_report(self) -> None:
        """리포트를 콘솔에 출력"""
        print(self.get_report())


def parse_arguments() -> argparse.Namespace:
    """
    명령줄 인자 파싱

    Returns:
        파싱된 인자 Namespace
    """
    parser = argparse.ArgumentParser(
        description="DTC53 Markers Parser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python parser.py --mode dtc_53 --sample-list samples.txt --plate-barcode PLATE_BARCODE
    python parser.py --mode irs_v2 --sample-list samples.txt --plate-barcode PLATE_BARCODE
        """
    )

    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level (default: INFO)"
    )

    parser.add_argument(
        "--conf-json",
        type=str,
        default="conf.json",
        help="Settings JSON path (default: conf.json)"
    )

    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        help="Mode key in settings JSON"
    )

    parser.add_argument(
        "--sample-list",
        type=str,
        required=True,
        default=None,
        help="Sample list file path (override config)"
    )

    parser.add_argument(
        "--plate-barcode",
        type=str,
        required=True,
        default="PLATE",
        help="Plate barcode to write into output OA file"
    )

    return parser.parse_args()


def main() -> int:
    """
    메인 함수

    Returns:
        종료 코드 (0: 성공, 1: 실패)
    """
    args = parse_arguments()

    try:
        conf = load_conf(args.conf_json, args.mode)
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}", file=sys.stderr)
        return 1

    apt_vcf = conf["apt_vcf"]
    imputed_vcf = conf["imputed_vcf"]
    rsmarker_list = conf["rsmarker_list"]
    sample_list_path = args.sample_list  # argparse 필수 입력
    plate_barcode = args.plate_barcode
    output_dir = Path(conf["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    output_oa_file = output_dir / f"{args.mode}.oa.txt"
    qc_report_path = conf.get("qc_report") or f"{output_oa_file}.qc.report.txt"
    Path(qc_report_path).parent.mkdir(parents=True, exist_ok=True)

    logger = setup_logger("parser", args.log_level, log_file=str(qc_report_path))
    logger.info("=" * 60)
    logger.info(f"{args.mode} Markers Parser Started")
    logger.info("=" * 60)

    try:
        qc = QCReporter(logger)
        qc.add_section(f"{args.mode} Markers Parser QC Report")
        qc.add_line(f"  실행 시간: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")

        # 1. 샘플 리스트 로드 및 분석
        sample_list_raw = load_list_file(sample_list_path, logger)
        sample_analysis = qc.analyze_list_file(sample_list_raw, "Sample List")
        sample_list = sample_analysis["unique_items"]  # 중복 제거된 리스트 사용

        # 2. RS 마커 리스트 로드 (TSV 형식) 및 분석
        marker_df = load_marker_tsv(rsmarker_list, logger)
        marker_list = marker_df['RS_ID'].tolist()
        qc.analyze_list_file(marker_list, "RS Marker List")

        # 마커 TSV 정보 출력
        qc.add_section("RS Marker TSV 정보")
        qc.add_key_value("총 마커 수", len(marker_df))
        qc.add_line("")
        qc.add_line("  [마커별 Allele 정보 (처음 10개)]")
        for idx, row in marker_df.head(10).iterrows():
            qc.add_line(f"    - {row['RS_ID']}: NORMAL={row['NORMAL_ALLELE']}, RISK={row['RISK_ALLELE']}")
        if len(marker_df) > 10:
            qc.add_line(f"    ... 외 {len(marker_df) - 10}개")

        # 3. APT VCF 파싱
        qc.add_section("APT VCF 파싱")
        qc.add_key_value("파일", apt_vcf)

        apt_parser = VCFParser(apt_vcf, logger)
        apt_parser.parse()

        apt_samples = apt_parser.get_samples()
        qc.add_key_value("VCF 내 샘플 수", len(apt_samples))
        qc.add_key_value("VCF 내 variant 수", len(apt_parser.data))

        # 4. APT VCF에서 마커 추출
        qc.add_section("APT VCF에서 마커 추출")
        apt_results, marker_positions = apt_parser.extract_markers(marker_df, sample_list)

        qc.add_key_value("요청 마커 수", len(marker_list))
        qc.add_key_value("APT VCF에서 찾은 마커 수", len(apt_results))

        # 찾지 못한 마커 (전체 출력)
        missing_in_apt = [m for m in marker_list if m not in apt_results]
        if missing_in_apt:
            qc.add_line("")
            qc.add_line(f"  [APT VCF에 없는 마커: {len(missing_in_apt)}개]")
            for m in missing_in_apt:
                qc.add_line(f"    - {m}")
            logger.info(f"APT VCF에 없는 마커 목록({len(missing_in_apt)}개): {', '.join(missing_in_apt)}")

        # 5. Imputed VCF 파싱 (헤더만) 및 필요 위치 수집
        qc.add_section("Imputed VCF 파싱")
        qc.add_key_value("파일", imputed_vcf)

        imputed_parser = VCFParser(imputed_vcf, logger)
        imputed_parser.parse_header_only()

        imputed_samples = imputed_parser.get_samples()
        qc.add_key_value("VCF 내 샘플 수", len(imputed_samples))

        needed_positions = set()
        for _, samples_data in apt_results.items():
            for _, geno_info in samples_data.items():
                if geno_info['is_no_call']:
                    chrom = str(geno_info['chrom']).replace('chr', '')
                    pos = str(geno_info['pos'])
                    needed_positions.add((chrom, pos))

        # APT에 아예 없는 RS들도 imputed에서 찾기 위해 position 수집 (RS는 imputed에 없음)
        missing_positions = {}
        if {'CHROM', 'POS'}.issubset(set(marker_df.columns)):
            for _, row in marker_df[marker_df['RS_ID'].isin(missing_in_apt)].iterrows():
                chrom = str(row['CHROM']).replace('chr', '')
                pos = str(row['POS'])
                missing_positions[row['RS_ID']] = (chrom, pos)

        # Imputed에서 조회할 위치 = APT no-call 위치 + (APT에 없는 RS의 위치)
        missing_position_set = set(missing_positions.values())
        total_positions = needed_positions | missing_position_set

        qc.add_key_value("Imputed 조회 위치 수", len(total_positions))

        if total_positions:
            qc.add_line("")
            qc.add_line("  [Imputed 조회 위치 목록 (chrom:pos)]")
            for chrom, pos in sorted(total_positions, key=lambda x: (str(x[0]), int(x[1]) if str(x[1]).isdigit() else str(x[1]))):
                qc.add_line(f"    - {chrom}:{pos}")
            logger.info(
                "Imputed 조회 위치 목록: "
                + ", ".join([f"{c}:{p}" for c, p in sorted(total_positions, key=lambda x: (str(x[0]), str(x[1])))])
            )

        # Imputed VCF에서 필요한 위치/RS만 인덱싱
        imputed_parser.index_positions(
            total_positions,
            target_samples=set(sample_list),
            missing_rs_ids=set()  # imputed에 RS ID가 없으므로 RS 기반 조회는 비활성화
        )
        qc.add_key_value("조회된 variant 수", imputed_parser.indexed_variant_count)

        # 6. No call 채우기 (imputed VCF에서)
        qc.add_section("No Call 처리")
        no_call_count = 0
        filled_from_imputed = 0
        filled_from_normal = 0

        for rs_id, samples_data in apt_results.items():
            for sample_id, geno_info in samples_data.items():
                if geno_info['is_no_call']:
                    no_call_count += 1

                    # imputed VCF에서 같은 position 찾기
                    chrom = geno_info['chrom']
                    pos = geno_info['pos']

                    imputed_variant = imputed_parser.get_variant_by_position(chrom, pos)

                    if imputed_variant is not None and sample_id in imputed_samples:
                        gt = imputed_variant[sample_id]
                        ref = imputed_variant['REF']
                        alt = imputed_variant['ALT']

                        a1, a2 = imputed_parser.decode_genotype(gt, ref, alt)

                        if a1 is not None and a2 is not None:
                            geno_info['allele1'] = a1
                            geno_info['allele2'] = a2
                            geno_info['is_no_call'] = False
                            geno_info['source'] = 'IMPUTED'
                            filled_from_imputed += 1
                            continue

                    # imputed에도 없으면 NORMAL_ALLELE로 채우기
                    marker_info = marker_df[marker_df['RS_ID'] == rs_id]
                    if len(marker_info) > 0:
                        normal_allele = marker_info.iloc[0]['NORMAL_ALLELE']
                        geno_info['allele1'] = normal_allele
                        geno_info['allele2'] = normal_allele
                        geno_info['is_no_call'] = False
                        geno_info['source'] = 'NORMAL_FILLED'
                        filled_from_normal += 1

        # 7. APT에 없는 마커도 imputed로 우선 채우기 (RS 대신 position 사용)
        imputed_filled_missing = 0
        for rs_id in missing_in_apt:
            pos_info = missing_positions.get(rs_id)
            imputed_variant = None
            if pos_info:
                chrom, pos = pos_info
                imputed_variant = imputed_parser.get_variant_by_position(chrom, pos)

            if imputed_variant is not None:
                chrom = imputed_variant['CHROM']
                pos = imputed_variant['POS']
                ref = imputed_variant['REF']
                alt = imputed_variant['ALT']

                apt_results[rs_id] = {}
                for sample_id in sample_list:
                    if sample_id in imputed_samples:
                        gt = imputed_variant.get(sample_id, './.')
                        a1, a2 = imputed_parser.decode_genotype(gt, ref, alt)
                        if a1 is not None and a2 is not None:
                            apt_results[rs_id][sample_id] = {
                                'allele1': a1,
                                'allele2': a2,
                                'chrom': chrom,
                                'pos': pos,
                                'ref': ref,
                                'alt': alt,
                                'gt_raw': gt,
                                'is_no_call': False,
                                'source': 'IMPUTED'
                            }
                            imputed_filled_missing += 1
                            continue

                    # imputed에 없거나 no-call이면 NORMAL_ALLELE로 채우기
                    marker_info = marker_df[marker_df['RS_ID'] == rs_id]
                    if len(marker_info) > 0:
                        normal_allele = marker_info.iloc[0]['NORMAL_ALLELE']
                        apt_results[rs_id][sample_id] = {
                            'allele1': normal_allele,
                            'allele2': normal_allele,
                            'chrom': '.',
                            'pos': '.',
                            'ref': normal_allele,
                            'alt': '.',
                            'gt_raw': 'FILLED',
                            'is_no_call': False,
                            'source': 'NORMAL_FILLED'
                        }
            else:
                # 위치 정보가 없거나 imputed에도 없으면 NORMAL_ALLELE로 채우기
                marker_info = marker_df[marker_df['RS_ID'] == rs_id]
                if len(marker_info) > 0:
                    normal_allele = marker_info.iloc[0]['NORMAL_ALLELE']
                    apt_results[rs_id] = {}
                    for sample_id in sample_list:
                        apt_results[rs_id][sample_id] = {
                            'allele1': normal_allele,
                            'allele2': normal_allele,
                            'chrom': '.',
                            'pos': '.',
                            'ref': normal_allele,
                            'alt': '.',
                            'gt_raw': 'FILLED',
                            'is_no_call': False,
                            'source': 'NORMAL_FILLED'
                        }

        qc.add_key_value("총 No Call 수", no_call_count)
        qc.add_key_value("Imputed VCF에서 채운 수", filled_from_imputed)
        qc.add_key_value("NORMAL_ALLELE로 채운 수", filled_from_normal)
        qc.add_key_value("Imputed로 채운 누락 RS 수", imputed_filled_missing)

        # 8. APT VCF에 없는 마커 중 여전히 없거나 no-call인 경우 NORMAL_ALLELE로 채우기
        for rs_id in missing_in_apt:
            if rs_id not in apt_results:
                marker_info = marker_df[marker_df['RS_ID'] == rs_id]
                if len(marker_info) > 0:
                    normal_allele = marker_info.iloc[0]['NORMAL_ALLELE']
                    apt_results[rs_id] = {}
                    for sample_id in sample_list:
                        apt_results[rs_id][sample_id] = {
                            'allele1': normal_allele,
                            'allele2': normal_allele,
                            'chrom': '.',
                            'pos': '.',
                            'ref': normal_allele,
                            'alt': '.',
                            'gt_raw': 'FILLED',
                            'is_no_call': False,
                            'source': 'NORMAL_FILLED'
                        }

        # 8. OA 형식으로 변환
        qc.add_section("OA 형식 출력")

        oa_rows = []
        for rs_id in marker_list:  # 원래 순서 유지
            if rs_id in apt_results:
                for sample_id in sample_list:
                    if sample_id in apt_results[rs_id]:
                        geno_info = apt_results[rs_id][sample_id]
                        source = geno_info.get('source', 'APT')

                        oa_row = {
                            'Sample ID': sample_id,
                            'Plate Barcode': plate_barcode,
                            'Gene Symbol': source,  # 이전 source를 Gene Symbol에 기록
                            'NCBI SNP Reference': rs_id,
                            'Assay Name or ID': rs_id,
                            'Allele 1 Call': geno_info['allele1'],
                            'Allele 2 Call': geno_info['allele2']
                        }
                        oa_rows.append(oa_row)

        result_df = pd.DataFrame(oa_rows)

        # 샘플명 우선 정렬 (샘플 리스트 순서 유지), 마커는 원래 순서 유지
        result_df['Sample ID'] = pd.Categorical(result_df['Sample ID'], categories=sample_list, ordered=True)
        result_df['NCBI SNP Reference'] = pd.Categorical(result_df['NCBI SNP Reference'], categories=marker_list, ordered=True)
        result_df = result_df.sort_values(['Sample ID', 'NCBI SNP Reference']).reset_index(drop=True)

        # JSON 형식 생성
        json_rows = []
        for rs_id in marker_list:
            if rs_id not in apt_results:
                continue
            for sample_id in sample_list:
                if sample_id not in apt_results[rs_id]:
                    continue
                geno_info = apt_results[rs_id][sample_id]
                a1 = sanitize_allele(geno_info['allele1'])
                a2 = sanitize_allele(geno_info['allele2'])
                json_rows.append({
                    "sampleId": sample_id,
                    "plateBarcode": plate_barcode,
                    "itemCd": rs_id,
                    "result1": a1,
                    "result2": a2
                })

        qc.add_key_value("총 출력 row 수", len(result_df))
        qc.add_key_value("샘플 수", len(sample_list))
        qc.add_key_value("마커 수", len(marker_list))

        # Source별 통계
        if 'Plate Barcode' in result_df.columns:
            source_counts = result_df['Plate Barcode'].value_counts().to_dict()
            qc.add_line("")
            qc.add_line("  [데이터 소스별 통계]")
            for source, count in source_counts.items():
                qc.add_line(f"    - {source}: {count} rows")

        # 9. OA 형식으로 파일 출력
        output_path = output_oa_file
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # JSON 파일 경로
        output_json_path = output_path.with_suffix(".json")

        # 메타데이터 생성
        meta = {
            'Source': f'{args.mode} Markers Parser',
            'Mode': args.mode,
            'APT VCF': apt_vcf,
            'Imputed VCF': imputed_vcf,
            'Marker List': rsmarker_list,
            'Sample List': sample_list_path,
            'Plate Barcode': plate_barcode,
            'Generated Date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
        }

        with open(output_path, 'w', encoding='utf-8') as f:
            for key, value in meta.items():
                f.write(f"# {key} : {value}\n")
            f.write("\n")
            result_df.to_csv(f, sep='\t', index=False)

        logger.info(f"Written {len(result_df)} rows to {output_path}")

        # JSON 파일 출력
        with open(output_json_path, 'w', encoding='utf-8') as jf:
            json.dump(json_rows, jf, ensure_ascii=False, indent=2)
        logger.info(f"Written {len(json_rows)} JSON records to {output_json_path}")

        # 10. QC 리포트 출력
        qc.add_section("처리 완료")
        qc.add_key_value("출력 파일", str(output_oa_file))
        qc.add_key_value("상태", "SUCCESS")

        # 콘솔에 리포트 출력
        qc.print_report()

        # 파일로 리포트 저장 (항상 저장, 경로 기본값 제공)
        qc_output_path = qc_report_path
        qc.write_report(qc_output_path)

        logger.info("=" * 60)
        logger.info("Process finished successfully")
        logger.info("=" * 60)
        return 0

    except Exception as e:
        logger.exception(f"Error during processing: {e}")
        logger.info("=" * 60)
        logger.info("Process finished with errors")
        logger.info("=" * 60)
        return 1


if __name__ == "__main__":
    sys.exit(main())

