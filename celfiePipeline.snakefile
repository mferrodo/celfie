configfile: "config.yaml"

# Sample list
baseIDS, = glob_wildcards("../bismarkcov/{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz")

rule all:
    input:
        "celfie_runFiles/testset.celfie",
        "celfieOut/1_alpha.pkl",
        "celfieOutput_deconvolution.tsv"

if config["build"] == "hg19":
    rule convert:
        input:
            cov = "../bismarkcov/{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz"
        output:
            "celfie_runFiles/{sample}_GRCh38_chr.bed"
        params:
            liftOver = config["liftover"],
            chain = config["chain"],
        shell:
            """
                mkdir -p tmp
                zcat {input} | awk '{{print $1 "\t" $2 "\t" $3+1 "\t" $5 "\t" $5+$6}}' - > tmp/{wildcards.sample}.bed.tmp ;
                {params.liftOver} tmp/{wildcards.sample}.bed.tmp {params.chain} tmp/{wildcards.sample}_GRCh38_chr.bed.tmp  tmp/unmapped.tsv.tmp ;
                sed 's/^/chr/' tmp/{wildcards.sample}_GRCh38_chr.bed.tmp  > {output} ;
                rm tmp/{wildcards.sample}.bed.tmp tmp/{wildcards.sample}_GRCh38_chr.bed.tmp
            """
else:
    rule convert:
        input:
            cov = "../bismarkcov/{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz"
        output:
            "celfie_runFiles/{sample}_GRCh38_chr.bed"
        shell:
            """
				zcat {input} | awk '{{print $1 "\t" $2 "\t" $3+1 "\t" $5 "\t" $5+$6}}' -  > {output}
            """

rule prepareCelfie:
    input:
        "celfie_runFiles/{sample}_GRCh38_chr.bed"
    output:
        "celfie_runFiles/{sample}_tims.txt"
    params:
        sites=config["sites"],
        sumbylist=config["sumbylist"],
    shell:
        """
        mkdir -p tmp;
		bedtools intersect -a "{input}" -b "{params.sites}" > tmp/{wildcards.sample}_500.txt.tmp ;
        python3 "{params.sumbylist}" "{params.sites}" tmp/{wildcards.sample}_500.txt.tmp {output} 1
        """

rule mergeCelfie:
    input:
        expand("celfie_runFiles/{sample}_tims.txt", sample = baseIDS)
    output:
        tims = "celfie_runFiles/testset.celfie",
        sample_order = "celfie_runFiles/testset_order.celfie"
    params:
        reference_tims_txt = config["reference_tims_summed"]
    run:
        import pandas as pd
        import numpy as np
        import sys
        import os
        import glob

        orderlist = []
        df_merged = pd.DataFrame()

        for file in input:
                print(file)
                file_name = os.path.splitext(os.path.basename(file))[0]
                df = pd.read_csv(file, sep="\t", header=None)
                orderlist.append(file_name)
                if df_merged.empty == True:
                    df_merged = df
                elif df_merged.empty == False:
                    df_merged = pd.merge(df_merged, df, how = "inner", right_on = [0,1,2], left_on = [0,1,2])

        reference_tims = pd.read_csv(params.reference_tims_txt, sep = "\t", header = None)

        df_final =  pd.concat([df_merged, reference_tims], axis = 1)

        df_final.to_csv(output.tims, sep = "\t", header = None, index = None)
        pd.DataFrame(orderlist).transpose().to_csv(output.sample_order, sep = "\t", header = None, index = None)

rule runCelfie:
    input:
        tims_testset = "celfie_runFiles/testset.celfie",
        sample_order = "celfie_runFiles/testset_order.celfie"
    output:
        pkl = "celfieOut/1_alpha.pkl",
        outdir = directory("celfieOut")
    params:
        celfie_em = config["celfie_em"],
        num_samples = len(baseIDS),
        iterations = config["iterations"],
        num_unk = config["num_unk"],
        iteration_number = config["iteration_number"],
        convergence_criteria = config["convergence_criteria"], 
        num_random_restart = config["num_random_restart"],
        random_seed = config["random_seed"]
    shell:
        """
        echo "Number of samples:"
        echo {params.num_samples}
        python3 {params.celfie_em} "{input.tims_testset}" "{output.outdir}" {params.num_samples} {params.iterations} {params.num_unk} {params.iteration_number} {params.convergence_criteria} {params.num_random_restart} {params.random_seed}
        """


rule labelCelfie:
    input:
        pkl = "celfieOut/1_alpha.pkl",
        sample_order = "celfie_runFiles/testset_order.celfie"
    output:
        celfieOutputFile = "celfieOutput_deconvolution.tsv"
    params:
        referenceLabels = config['referenceLabels']
    run:
        import pickle
        import csv
        import os
        import pandas as pd

        infile = open(input.pkl,'rb')
        dict = pickle.load(infile)
        df = pd.DataFrame.from_dict(dict)

        inputRefKey = params.referenceLabels
        with open(inputRefKey, 'r') as f:
          reader = csv.reader(f, delimiter="\t")
          column_list = list(reader)
        flat_list = [item for sublist in column_list for item in sublist]
        df.columns = flat_list

        inputSampleKey = input.sample_order
        with open(inputSampleKey, 'r') as f:
          reader = csv.reader(f, delimiter="\t")
          column_list = list(reader)

        flat_list = [item for sublist in column_list for item in sublist]

        df.index = flat_list
        deconvOut = output.celfieOutputFile
        df.to_csv(deconvOut, sep = "\t")
