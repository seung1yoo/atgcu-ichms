python -m atgcu_ichms.cli extract \
    --product irs_v2 \
    --path-anl-for-all data/251204_IPMI_Analysis_For_All/ \
    --path-wkdir data/251204_IPMI_Analysis_For_irs_v2 \
    --sample-sheet data/251204_IPMI.sample_sheet.tsv \
    --imputed-mode tabix

python -m atgcu_ichms.cli api-to-lis \
    --product irs_v2 \
    --json data/251204_IPMI_Analysis_For_irs_v2/Analysis_For_irs_v2/irs_v2.oa.json \
    --insecure
