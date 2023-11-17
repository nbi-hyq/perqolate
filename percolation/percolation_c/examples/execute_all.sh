for exe in $(file * | grep executable | cut -f1 -d:)
do
  nohup stdbuf -oL ./$exe > $exe.txt &
done
wait
