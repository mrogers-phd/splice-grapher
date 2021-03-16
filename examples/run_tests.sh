#! /bin/bash
# Copyright (C) 2010 by Colorado State University
# Contact: Mark Rogers <rogersma@cs.colostate.edu>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
# USA.
function checkErrs {
    if [ $1 -eq 1 ] ; then
        echo "** Error in $2"
        exit 1
    fi
}

function logmsg {
    echo $1
    echo $1 >> run_tests.log
}

allGenes=""
logFile=run_tests.log
logmsg "== Converting gene models =="
for f in $(ls *_model.gff); do
    gene=$(echo $f | sed 's/_model\.gff//')
    logmsg "  converting $f to splicegraph"
    logmsg "    gene_model_to_splicegraph.py -m $f -g $gene -a -o ${gene}_graph.gff"
    gene_model_to_splicegraph.py -m $f -g $gene -a -o ${gene}_graph.gff
    checkErrs $? gene_model_to_splicegraph.py
    allGenes="$gene ${allGenes}"
done

logmsg ""
logmsg "== Converting EST alignments =="
for f in $(ls *_ests.psl); do
    gene=$(echo $f | sed 's/_ests\.psl//')
    model=${gene}_model.gff
    logmsg "  converting $f to splicegraph"
    logmsg "    ests_to_splicegraph.py $f -i 10 -g $gene -m ${model}"
    ests_to_splicegraph.py $f -i 10 -g $gene -m ${model}
    checkErrs $? ests_to_splicegraph.py
done

logmsg ""
logmsg "== Generating predictions =="
outputPdf=""
for g in ${allGenes}; do
    model=${g}_model.gff
    modelgraph=${g}_graph.gff
    ests=${g}_ests.gff
    reads=${g}.sam
    prediction=${g}_predicted.gff
    plotfile=${g}_predicted.pdf
    cfgFile=${g}_plot.cfg
    plotterFile=${g}_plot.pdf

    for f in ${modelgraph} ${ests} ${reads} ${known}; do
        if [ ! -e $f ]; then
            logmsg "** Error: cannot locate $f; skipping $g"
            continue
        fi
    done

    pred_params="-d ${reads} -s ${ests} -o ${prediction} -J 2 -M 10 -T 2"
    logmsg "  predict_splicegraph.py ${modelgraph} ${pred_params}"
    predict_splicegraph.py ${modelgraph} ${pred_params}
    checkErrs $? predict_splicegraph.py

    if [ ! -e $prediction ]; then
        logmsg "** Error: no prediction produced for ${g}"
        continue
    fi

    ##view_params="-d ${reads} -m ${model} -o ${plotfile} -s ${prediction} -G ${modelgraph} -J 2 -cL"
    ##logmsg "  view_splicegraph_multiplot.py ${g} ${view_params}"
    ##view_splicegraph_multiplot.py ${g} ${view_params}
    ##checkErrs $? view_splicegraph_multiplot.py
    ##if [ ! -e $plotfile ]; then
    ##    logmsg "** Error: no plot produced for ${g}"
    ##else
    ##    outputPdf="$plotfile $outputPdf"
    ##fi

    view_params="-d ${reads} -m ${model} -o ${plotfile} -s ${prediction} -G ${modelgraph} -J 2 -cL"
    logmsg "  plotter.py ${cfgFile}"
    plotter.py ${cfgFile}
    checkErrs $? plotter.py

    if [ ! -e $plotterFile ]; then
        logmsg "** Error: no plotter output produced for ${g}"
    else
        outputPdf="$plotterFile $outputPdf"
    fi
done

logmsg ""
if [[ "$outputPdf" != "" ]]; then
    logmsg "Generated the following plot files:"
    for f in $outputPdf; do
        logmsg "  $f"
    done
fi
