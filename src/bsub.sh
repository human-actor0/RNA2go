
bsub.view(){
OUTPUT=${@:$#};
cmd="
#!/bin/sh
#BSUB -n 1
#BSUB -o $OUTPUT.bsub.out.%J
#BSUB -e $OUTPUT.bsub.err.%J
#BSUB -J $OUTPUT
hm util mkdir $OUTPUT"
"
echo "$cmd";
echo "${@:-1}"
#echo "$cmd" >$OUTPUT.cmd;
#echo $@ >> $OUTPUT.cmd;
#bsub < $OUTPUT.cmd
}
