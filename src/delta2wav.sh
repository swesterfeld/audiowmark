for delta in 0.080 0.050 0.030 0.020 0.015 0.010 0.005
do
  audiowmark add "$1" --water-delta=$delta /tmp/w$delta.wav a9fcd54b25e7e863d72cd47c08af46e61
done

audiowmark scale "$1" /tmp/w.orig.wav
