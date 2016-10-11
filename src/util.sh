

util.split(){
usage="$FUNCNAME <file> <by_nlines> <out> 
"
if [ $# -lt 3 ];then echo "$usage"; return; fi
	awk -v n=$2 -v o=$3 'NR%n==1{x=sprintf("%s.%d",o,++i);}{print > x}' $1
#	split -l $2 $1 $3
}
util.split.test(){
echo \
"A	a
B	b
C	c
D d
E f	g" > tmp.i
util.split tmp.i 2 tmp.o
head tmp.o*
rm tmp.i tmp.o*

}

util.prefixsufix(){
	EXPANDED=()
	for E in `cat ${@:3}`; do
	    EXPANDED+=("$1${E}$2")
	done
	echo "${EXPANDED[@]}"
}
util.prefixsufix.test(){
	echo "A B C" | util.prefixsufix "pre" "suf" -
}
util.mkdir(){
	if [ -z "${1##*\/*}" ];then
                mkdir -p ${1%/*}
	fi
}
util.mkout(){
	local res=$1;
        if [ ! -z ${1##*${2}} ];then 
                res=$res$2; 
        fi   
	util.mkdir $res
	echo "$res"
}

util.mktempd(){
	mktemp -d 2>/dev/null || mktemp -d -t 'hmtmpdir'
}
util.check(){ ## obtained from bamtools test
	if diff $1 $2; then
		echo ok
		return 1
	else
		echo fail
		return 0
	fi
}
util.run_R(){
usage="
usage: $FUNCNAME <Rscript> <out> 
"
	local tmpd=`hm util mktempd`;
	local cmd=$@;
	cmd=${cmd/stdout/$tmpd\/out};
	echo "$cmd" > $tmpd/cmd;

        R --no-save -q -f $tmpd/cmd &> $tmpd/log;
	if [ -f $tmpd/out ];then
		cat $tmpd/out
	else
		cat $tmpd/log
	fi
	local flag=${2:-""};

	if [ "$flag" = "log" ];then 
		cat $tmpd/log
	fi
	rm -rf $tmpd;
}
