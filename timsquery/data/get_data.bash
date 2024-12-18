# LFQ_timsTOFPro_diaPASEF_Ecoli_03.d.zip

# Long gradient run
docker run --rm --platform linux/amd64 -v ${PWD}:/data \
  ghcr.io/pride-archive/aspera \
  ascp -i /home/aspera/.aspera/cli/etc/asperaweb_id_dsa.openssh \
  -TQ -P33001 \
  prd_ascp@fasp.ebi.ac.uk:/pride/data/archive/2022/02/PXD028735/LFQ_timsTOFPro_diaPASEF_Condition_A_Sample_Alpha_02.d.zip \
  /data

# 22 min hela from some random day
aws s3 cp s3://terraform-workstations-bucket/jspaezp/20241027_PRTC/230510_PRTC_13.d.tar .
tar -xf 230510_PRTC_13.d.tar
