
util.mktempd(){
	mktemp -d 2>/dev/null || mktemp -d -t 'hmtmpdir'
}
