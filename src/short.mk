all: short-60 short-60-mp3 \
  short-50 short-50-mp3 \
  short-40 short-40-mp3 \
  short-30 short-30-mp3 \
  short-20 short-20-mp3 \
  short-10 short-10-mp3

short-60:
	AWM_FILE=t-short-60 AWM_CLIP=60 fer-test.sh 10 "" > tmp-short-60
	mv tmp-short-60 short-60

short-60-mp3:
	AWM_FILE=t-short-60-mp3 AWM_CLIP=60 fer-test.sh 10 "" mp3 128 > tmp-short-60-mp3
	mv tmp-short-60-mp3 short-60-mp3

short-50:
	AWM_FILE=t-short-50 AWM_CLIP=50 fer-test.sh 10 "" > tmp-short-50
	mv tmp-short-50 short-50

short-50-mp3:
	AWM_FILE=t-short-50-mp3 AWM_CLIP=50 fer-test.sh 10 "" mp3 128 > tmp-short-50-mp3
	mv tmp-short-50-mp3 short-50-mp3

short-40:
	AWM_FILE=t-short-40 AWM_CLIP=40 fer-test.sh 10 "" > tmp-short-40
	mv tmp-short-40 short-40

short-40-mp3:
	AWM_FILE=t-short-40-mp3 AWM_CLIP=40 fer-test.sh 10 "" mp3 128 > tmp-short-40-mp3
	mv tmp-short-40-mp3 short-40-mp3

short-30:
	AWM_FILE=t-short-30 AWM_CLIP=30 fer-test.sh 10 "" > tmp-short-30
	mv tmp-short-30 short-30

short-30-mp3:
	AWM_FILE=t-short-30-mp3 AWM_CLIP=30 fer-test.sh 10 "" mp3 128 > tmp-short-30-mp3
	mv tmp-short-30-mp3 short-30-mp3

short-20:
	AWM_FILE=t-short-20 AWM_CLIP=20 fer-test.sh 10 "" > tmp-short-20
	mv tmp-short-20 short-20

short-20-mp3:
	AWM_FILE=t-short-20-mp3 AWM_CLIP=20 fer-test.sh 10 "" mp3 128 > tmp-short-20-mp3
	mv tmp-short-20-mp3 short-20-mp3

short-10:
	AWM_FILE=t-short-10 AWM_CLIP=10 fer-test.sh 10 "" > tmp-short-10
	mv tmp-short-10 short-10

short-10-mp3:
	AWM_FILE=t-short-10-mp3 AWM_CLIP=10 fer-test.sh 10 "" mp3 128 > tmp-short-10-mp3
	mv tmp-short-10-mp3 short-10-mp3
