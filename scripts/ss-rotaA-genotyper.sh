
deploy_libdir=$KB_TOP/lib/rota-a-genotyper
dev_libdir=$KB_TOP/modules/bvbrc_subspecies_classification/lib/rota-a-genotyper

if [[ -d $deploy_libdir ]] ; then
    libdir=$deploy_libdir
else
    if [[ -d $dev_libdir ]] ; then
	libdir=$dev_libdir
    else
	echo "Cannot find genotyper libdir in $deploy_libdir or $dev_libdir" 1>&2
	exit 1
    fi
fi

java -jar $libdir/StandAloneRtvAGenotyper.jar $libdir/rotaAGenotyper.config "$@"
       
    
