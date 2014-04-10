#!/usr/bin/env bash
LIBDIR=libseqlink
CMDS=seqlink
PROGRAM=SEQLinkage
set -e
is_relative() {
    local path="$1"
    shift
    [ "${path:0:1}" != "/" ]
    return
}
install () {
    mkdir -p $2/bin
    rm -rf $2/lib/$LIBDIR &> /dev/null
    mkdir -p $2/lib/$LIBDIR
    cp -a $1/* $2/lib/$LIBDIR
    for cmd in $CMDS; do
	    rm -f $2/bin/$cmd &> /dev/null
        cd $2/bin
	    ln -s ../lib/$LIBDIR/$cmd $cmd
    done
    echo -e "Libraries are installed to $2/lib\nBinary files are installed to $2/bin\n"
}
main () {
    local fullpath=$2
    if [ -z $fullpath ]; then
        echo "Enter installation directory for $PROGRAM: "
        printf "\t [/usr/local]  "
        read fullpath
    fi
    if [ -z $fullpath ]; then
	    install $1 "/usr/local"
    else
	    eval fullpath=$fullpath
	    if is_relative $fullpath; then
	        fullpath=$PWD/$fullpath
	    fi
	    install $1 $fullpath
    fi
}
main $@
