import pandas as pd
from sklearn.datasets import load_wine
from sklearn.linear_model import LogisticRegression

import flytekit.extras.sklearn
from flytekit import task, workflow
import os
from typing import Tuple

import flytekit
from flytekit import kwtypes, task, workflow
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile

t1 = ShellTask(
    name="task1",
    debug=True,
    script="""
    rm -f ali.sam
    ./minimap2 -a -t 6 GCF_000005845.2_ASM584v2_genomic.fna SRR24658890.fasta > {outputs.i}
    echo minimap2 finished
    """,
    output_locs=[
        OutputLocation(var="i", var_type=FlyteFile, location="ali.sam")
    ],
)

t2 = ShellTask(
    name="task2",
    debug=True,
    script="""
    rm -f {outputs.j}
    touch {outputs.j}
    QUALITY=$(samtools flagstat {inputs.x} | python3 -c 'from sys import stdin; f=stdin.read(); import re; print(re.findall("\d+.\d+%", f)[0].replace("%", ""))')
    echo $QUALITY >> {outputs.j}
    """,
    inputs=kwtypes(x=FlyteFile),
    output_locs=[
        OutputLocation(var="j", var_type=FlyteFile, location="quality.txt")
    ],
)

t3 = ShellTask(
    name="task2",
    debug=True,
    script="""
    python3 -c 'f = open("{inputs.x}", "r", encoding="utf-8").read();\nif int(float(f)) > 90: print("OK");\nelse: print("NOT OK")'
    rm -f {inputs.x}
    """,
    inputs=kwtypes(x=FlyteFile),
)

t4 = ShellTask(
    name="task2",
    debug=True,
    script="""
    rm -f {outputs.k}
    samtools view -S -b {inputs.x} > {outputs.k}
    """,
    inputs=kwtypes(x=FlyteFile),
    output_locs=[
        OutputLocation(var="k", var_type=FlyteFile, location="ali.bam")
    ],
)

t5 = ShellTask(
    name="task2",
    debug=True,
    script="""
    rm -f {outputs.l}
    samtools sort {inputs.x} -o {outputs.l}
    """,
    inputs=kwtypes(x=FlyteFile),
    output_locs=[
        OutputLocation(var="l", var_type=FlyteFile, location="ali.sorted.bam")
    ],
)

t6 = ShellTask(
    name="task2",
    debug=True,
    script="""
    rm -f {outputs.m}
    freebayes -f GCF_000005845.2_ASM584v2_genomic.fna {inputs.x} > {outputs.m}
    """,
    inputs=kwtypes(x=FlyteFile),
    output_locs=[
        OutputLocation(var="m", var_type=FlyteFile, location="var.vcf")
    ],
)


@workflow
def wf():
    ali = t1()
    q = t2(x=ali)
    t3(x=q)
    bam = t4(x=ali)
    ali_sorted = t5(x=bam)
    t6(x=ali_sorted)

if __name__ == "__main__":
    wf()

# rm -f ali.sam
#
# ./minimap2 -a -t 6 GCF_000005845.2_ASM584v2_genomic.fna SRR24658890.fasta > ali.sam
# echo minimap2
#
# QUALITY=$(samtools flagstat ali.sam | python3 -c 'from sys import stdin; f=stdin.read(); import re; print(re.findall("\d+.\d+%", f)[0].replace("%", ""))')
# echo QUALITY=$QUALITY
#
# CMP=$(python3 -c "print(int(float($QUALITY) > 90))")
# echo CMP=$CMP
#
# rm -f ali.bam
# samtools view -S -b ali.sam > ali.bam
#
# rm -f ali.sorted.bam
# samtools sort ali.bam  -o ali.sorted.bam
#
# freebayes -f GCF_000005845.2_ASM584v2_genomic.fna ali.sorted.bam > var.vcf
