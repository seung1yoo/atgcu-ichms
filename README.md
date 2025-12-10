# atgcu-ichms

Markers 추출 및 LIS API 전송 파이프라인 CLI.

## 주요 기능
- `extract`: APT/Imputed VCF에서 마커를 추출해 OA(txt)와 JSON을 생성.
- No-call을 Imputed VCF로 보완, 없으면 NORMAL_ALLELE로 채움.
- 샘플/마커 순서를 유지해 정렬된 OA 출력.
- allele 정규화(`N`, `*`, `-`, `None` → `.`).
- `api-to-lis`: 생성된 JSON을 conf.json에 지정된 API로 POST(타임아웃/재시도 지원).

## 요구사항
- Python 3.8+
- pandas, requests (pysam은 선택: tabix 인덱스 사용 시 성능 향상)

## 설정 예시 (`conf.json`)
```json
{
  "product": {
    "apt_vcf": "/path/GenotypeCall.vcf",
    "imputed_vcf": "/path/imputed_all.vcf.gz",
    "rsmarker_list": "/path/dtc_53_rs_markers.txt",
    "output_dir": "/path/output",
    "api": {
      "endpoint": "https://example.com/lis/ingest",
      "headers": {"Content-Type": "application/json"},
      "auth_token": "YOUR_TOKEN"
    }
  }
}
```

## 사용법
진입점: `parser.py` (또는 `python -m atgcu_ichms.cli ...`)

### 1) 추출 (OA + JSON 생성)
```bash
python parser.py extract \
  --product product \
  --sample-list /path/samples.txt \
  --plate-barcode PLATE123 \
  --conf-json conf.json \
  --log-level INFO
```
- 출력: `output_dir/product.oa.txt`, `output_dir/product.json`, QC 리포트 `product.oa.txt.qc.report.txt`

### 2) JSON을 LIS API로 전송
```bash
python parser.py api-to-lis \
  --product product \
  --json /path/output/product.json \
  --conf-json conf.json \
  --timeout 10 \
  --retries 3
```
- `conf.json`의 `product.api.endpoint`, `headers`, `auth_token`을 사용.

## 참고
- Imputed VCF는 bgzip + tabix 인덱스가 있으면 성능이 향상됩니다.
- `pysam` 미설치 시 자동으로 스트리밍 fallback을 사용합니다.

