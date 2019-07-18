#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

// parameters
params.exon = ""
params.ref  = ""
params.assembly = "GRCh38"
params.vepCache = "/sc1/groups/bfx-red/data/datainsights/users/yaol12/vep_data/"  
params.vepDocker = "ensemblorg/ensembl-vep:release_96.0"
params.vepCacheVersion = "96"
params.outDir = "./"
params.codeDir = "/sc1/groups/bfx-red/users/yaol12/software/ecTMB_refbuild/inst/nf/"
params.mutContext96="/sc1/groups/bfx-red/data/datainsights/users/yaol12/GSMuta_data/GSMuta/mutation_context_96.txt"
params.nline = 10
params.pick = false
params.vcf2maf="/sc1/groups/bfx-red/users/yaol12/software/vcf2maf/tmp/vcf2maf/vcf2maf_edit.pl"


println "=======INFO==================================================="
println "=======BEGIN=================================================="
println "exon bed file: $params.exon"
println "ref path: $params.ref"
println "assembly: $params.assembly"
println "vep Cache version: $params.vepCacheVersion"
println "vepCache: $params.vepCache"
println "vepDocker:$params.vepDocker"
println "outDir:   $params.outDir"
println "codeDir:  $params.codeDir"
println "mutContext96: $params.mutContext96"
println "nline: $params.nline"
println "pick: $params.pick"
println "vcf2maf: $params.vcf2maf"
println "=======END===================================================="
println ""
println ""

bed2vcf = params.codeDir + "bed2vcf.changes.pl"
processVEP = params.codeDir + "process.vep.output.pl"
processMAF = params.codeDir + "process.maf.output.pl"
createRdata = params.codeDir + "create.exome.Rdata.R"

///////////////// Split bed file //////////////////////////////
process splitbed{
	tag ("$params.exon")
	validExitStatus 0,42
    cpus 1
    clusterOptions { "-q all.q -l h_vmem=${1 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }

    input:
    	file(exon) from Channel.fromPath(params.exon)


	output:
		file ('split_*') into splitbedCh mode flatten

	script:
		"""
		split -l ${params.nline} ${exon} split_
		"""

}


// process vepvcf{
// 	tag ("${bedN}")
// 	validExitStatus 0,42
//     cpus 1
//     clusterOptions { "-q all.q -l h_vmem=${5 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }

//     input:
//     	file(bed) from splitbedCh

//     output:
//     	file( "${bedN}.ann.vcf" ) into (annVcfCh, annVcfCh1, annVcfCh2)

//     script:
//     	bedN = bed.getName()
//         if(params.pick){
//             extra="--pick"
//         }else{
//             extra = ""
//         }

//     	"""
//     	perl ${bed2vcf} -e ${bed} -f ${params.ref}

//     	docker run --rm -u=\$UID -v /sc1/:/sc1/ -v ${params.vepCache}:/opt/vep/.vep  \
//     	-v \$PWD:\$PWD -w \$PWD ${params.vepDocker} /opt/vep/src/ensembl-vep/vep \
//     	--offline --no_stats --domains --failed 1 --total_length --biotype \
//     	--fork 1 --force_overwrite --fa ${params.ref} --dir /opt/vep/.vep \
//     	--vcf -i output.vcf -o ${bedN}.ann.vcf --assembly ${params.assembly} --cache_version ${params.vepCacheVersion} ${extra}
//     	"""
// }

annVcfCh1.toSortedList().set{annVcfChmerge}

process mergeVCF{
    tag("mergeVCF")
    validExitStatus 0,42
    cpus 1
    clusterOptions { "-q all.q -l h_vmem=${5 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }
    publishDir "${params.outDir}/${params.vepCacheVersion}_${params.assembly}/", mode: 'link'

    input:
        file "split_*.ann.vcf" from annVcfChmerge
    output:
        file "merged.ann.vcf" into finalVCF

    script:
        """
        ## output header
        grep "#" split_1.ann.vcf > header.tsv

        ## concat
        { cat header.tsv;  cat split_*.ann.vcf | grep -v '#'; } > merged.ann.vcf
        """
}



if(params.pick){
    process postprocess{
        tag("$vcfN")
        validExitStatus 0,42
        cpus 1
        clusterOptions { "-q all.q -l h_vmem=${5 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }

        input:
            file(vcf) from annVcfCh

        output:
            file( "${vcfN}_processed_vep.tsv" ) into annVEPCh

        script:
            vcfN = vcf.getName().replaceAll(/\.ann\.vcf/, "")

            """
            ## process vep, output name processed_vep.tsv
            perl ${processVEP} -i ${vcf}
            mv processed_vep.tsv ${vcfN}_processed_vep.tsv
            """
    }
    annVEPCh.toSortedList().set{annVEPChmerge}


    process mergeVEP{
        tag("mergeVEP")
        validExitStatus 0,42
        cpus 1
        clusterOptions { "-q all.q -l h_vmem=${5 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }
        publishDir "${params.outDir}/${params.vepCacheVersion}_${params.assembly}/", mode: 'link'

        input:
            file "split_*._processed_vep.tsv" from annVEPChmerge
        output:
            file "merged.processed_vep.tsv" into finalVEP

        script:
            """
            ## concat
            cat  split_*._processed_vep.tsv | sort -k 1,1 -k2,2n > merged.processed_vep.tsv
            """

    }

    process finalprocess{
        tag("merge")
        validExitStatus 0,42
        cpus 1
        clusterOptions { "-q highmem.q -l hp -l h_vmem=${190 + 50 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }
        publishDir "${params.outDir}/${params.vepCacheVersion}_${params.assembly}/", mode: 'link'

        input:
            file(vep) from finalVEP

        output:
            file("exome_${params.vepCacheVersion}_${params.assembly}_vep.Rdata") into finalch

        script:
            """
            ## 
            Rscript ${createRdata} ${vep} ${params.mutContext96} exome_${params.vepCacheVersion}_${params.assembly}_vep.Rdata
            """
    }

}else{
    
    process vcf_maf{
        tag("$vcfN")
        validExitStatus 0,42
        cpus 1
        clusterOptions { "-q all.q -l h_vmem=${5 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }

        input:
            file(vcf) from annVcfCh2

        output:
            file( "${vcfN}_processed_maf.tsv" ) into processedMAFCh

        script:
            vcfN = vcf.getName().replaceAll(/\.ann\.vcf/, "")

            """
            ## convert vcf to maf, 
            perl ${params.vcf2maf} --input-vcf ${vcf} --output-maf ${vcfN}.maf --ref-fasta  ${params.ref} --filter-vcf 0 --ncbi-build ${params.assembly} 
            paste <(grep -v '#' ${vcf} | cut -f8 | cut -d ';' -f 1 ) <( egrep -v "Hugo_Symbol|#"  ${vcfN}.maf) >  ${vcfN}.forprocess.maf
            perl ${processMAF} -i ${vcfN}.forprocess.maf
            mv processed_MAF.tsv ${vcfN}_processed_maf.tsv

            """

    }

    processedMAFCh.toSortedList().set{annMAFChmerge}

    process mergeMAF{
        tag("mergeMAF")
        validExitStatus 0,42
        cpus 1
        clusterOptions { "-q all.q -l h_vmem=${5 + 2 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }
        publishDir "${params.outDir}/${params.vepCacheVersion}_${params.assembly}/", mode: 'link'

        input:
            file "split_*._processed_maf.tsv" from annMAFChmerge
        output:
            file "merged.processed_maf.tsv" into finalMAF

        script:
            """
            ## concat
            cat  split_*._processed_maf.tsv | sort -k 1,1 -k2,2n > merged.processed_maf.tsv
            """

    }



    process finalprocessMAF{
        tag("merge")
        validExitStatus 0,42
        cpus 1
        clusterOptions { "-q highmem.q -l hp -l h_vmem=${190 + 50 * (task.attempt -1)}G -pe smp 1 -P rssrbfx -v PATH" }
        publishDir "${params.outDir}/${params.vepCacheVersion}_${params.assembly}/", mode: 'link'

        input:
            file(vep) from finalMAF

        output:
            file("exome_${params.vepCacheVersion}_${params.assembly}_vep.Rdata") into finalchmaf

        script:
            """
            ## 
            Rscript ${createRdata} ${vep} ${params.mutContext96} exome_${params.vepCacheVersion}_${params.assembly}_vep.Rdata
            """
    }

}











