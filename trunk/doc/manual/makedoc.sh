if [ ! $MARLINWORKDIR ] ; then
    export MARLINWORKDIR=$MARLIN
fi


latex manual.tex
latex manual.tex

dvips manual.dvi

latex2html -mkdir -dir $MARLINWORKDIR/doc/manual_html -info 0 -no_auto_link -split 0 -no_navigation manual.tex


