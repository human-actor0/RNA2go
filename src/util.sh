

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
                mkdir -p ${3%/*}
	fi
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
