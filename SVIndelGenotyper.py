#!/usr/bin/env python3

import os
import sys
import re
import subprocess
import json
import math
import random
import string

# Define constant paths and folders
blast_exe_folder="./tool/"


def print_usage():
    print(f"Usage: python {sys.argv[0]} [SV_list_fn] [optional: clean_deepseq_folder]")
    sys.exit()


def get_bit_fn(bam, reference_folder="./reference/"):
    bit_fn = reference_folder + "hg19.fa"
    genome_version = "hg19"
    chr_in_bam = "yes"
    genome_version_command = f"samtools idxstats {bam} | head -1"
    chr_checkresult = str(subprocess.check_output(genome_version_command, shell=True), "utf-8")
    if "chr1\t248956422" in chr_checkresult or "1\t248956422" in chr_checkresult:
        bit_fn = reference_folder + "hg38.fa"
        genome_version = "hg38"
    if chr_checkresult.count("chr") == 0:
        chr_in_bam = "no"
    return bit_fn, genome_version, chr_in_bam

def get_revcmp(strs):
    cpm = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "": "", "-": "-"}
    rev = strs[::-1]
    revcmp = "".join([cpm[i] for i in rev])
    return revcmp

def get_depth_indel(indel, sample, chr_in_bam, bam, clean_deepseq_folder="no"):
    chro, pos, ref_allele, mut_allele = indel.split(".")
    if clean_deepseq_folder != "no":
        clean_deepseq_fn = f"{clean_deepseq_folder}/{sample}.count.gz"
        depth_pos_command = f"tabix {clean_deepseq_fn} {chro.replace('chr','')}:{pos}-{pos}"
        depth_pos = "0"
        try:
            pos_result = str(subprocess.check_output(depth_pos_command, shell=True), "utf-8").strip()
            if len(pos_result) > 0:
                depth_pos = pos_result.split("\n")[0].split("\t")[-1]
        except:
            print(f"no depth: {chro}.{pos}")
    else:
        if chr_in_bam == "yes":
            depth_pos_command = f"samtools depth -r {chro}:{pos}-{pos} {bam} | cut -f3"
        else:
            depth_pos_command = f"samtools depth -r {chro.replace('chr', '')}:{pos}-{pos} {bam} | cut -f3"
        try:
            depth_pos = str(subprocess.check_output(depth_pos_command, shell=True), "utf-8").strip()
        except:
            print(f"no depth: {chro}.{pos}")   
    return depth_pos

def make_ref_indel(indel, bam, bit_fn):
    chro, pos, ref_allele, mut_allele = indel.split(".")
    sample = os.path.basename(bam)
    
    if "chr" not in indel:
        chro = "chr" + chro
    
    startA = int(pos) - 200
    endA = int(pos) - 1
    startB = int(pos) + len(ref_allele)
    endB = int(startB) + 200 - 1
    
    left_ref_command = f"samtools faidx {bit_fn} {chro}:{startA}-{endA}"
    right_ref_command = f"samtools faidx {bit_fn} {chro}:{startB}-{endB}"
    
    left_ref_seq = ("").join(str(subprocess.check_output(left_ref_command, shell=True), "utf-8").strip().split("\n")[1:]).upper()
    right_ref_seq = ("").join(str(subprocess.check_output(right_ref_command, shell=True), "utf-8").strip().split("\n")[1:]).upper()
    
    WT_ref_seq = left_ref_seq + ref_allele + right_ref_seq
    mut_ref_seq = left_ref_seq + mut_allele + right_ref_seq
    
    WT_genome_pos = [f"{chro}.{i}" for i in range(startA, endB + 1)]
    mut_genome_pos = [f"{chro}.{i}" for i in range(startA, endA + 2)] + [f"mut.{i}" for i in range(len(mut_allele) - 1)] + [f"{chro}.{i}" for i in range(startB, endB + 1)]
    
    WT_ref_fn = ''.join(random.choice(string.ascii_lowercase) for _ in range(7)) + "." + sample + "." + indel + "_WT_ref.fa"
    mut_ref_fn = ''.join(random.choice(string.ascii_lowercase) for _ in range(7)) + "." + sample + "." + indel + "_mut_ref.fa"
    
    with open(WT_ref_fn, "w") as WT_ref_fh:
        WT_ref_fh.write(">" + sample + "." + indel + "\n")
        WT_ref_fh.write(WT_ref_seq + "\n")
    
    with open(mut_ref_fn, "w") as mut_ref_fh:
        mut_ref_fh.write(">" + sample + "." + indel + "\n")
        mut_ref_fh.write(mut_ref_seq + "\n")
    
    ref_data = {
        "WT_ref_fn": WT_ref_fn,
        "mut_ref_fn": mut_ref_fn,
        "WT_ref_seq": WT_ref_seq,
        "mut_ref_seq": mut_ref_seq,
        "WT_genome_pos": WT_genome_pos,
        "mut_genome_pos": mut_genome_pos,
        "WT_keyseq": WT_ref_seq[(200 - 6):(200 + 6)],
        "WT_keyseq_pos": WT_genome_pos[(200 - 6):(200 + 6)],
        "mut_keyseq": mut_ref_seq[(200 - 6):(200 + len(mut_allele) + 6)],
        "mut_keyseq_pos": mut_genome_pos[(200 - 6):(200 + len(mut_allele) + 6)]
    }
    
    return ref_data

def extract_sam_indel(indel, bam, chr_in_bam):
    chro, pos, ref_allele, mut_allele = indel.split(".")
    sample = os.path.basename(bam).replace(".bam", "")
    get_chro_command = f"samtools idxstats {bam} | cut -f1"
    chr_checkresult = str(subprocess.check_output(get_chro_command, shell=True), "utf-8")
    
    if chr_in_bam == "yes" and "chr" not in indel:
        chro = "chr" + chro
    
    if chr_in_bam == "no" and "chr" in indel:
        chro = chro.replace("chr", "")
    
    sam_fn = ''.join(random.choice(string.ascii_lowercase) for _ in range(8)) + "." + indel + "." + os.path.basename(bam).replace(".bam", ".sam")
    command_pos = f"samtools view {bam} {chro}:{int(pos) - 100}-{int(pos) + 100} | sort -k1,1 > {sam_fn}"
    os.system(command_pos)
    
    return sam_fn

def get_read_qual(qual_str, qual_cutoff):
    read_len = len(qual_str)
    ngood = 0
    
    for nn_qual in qual_str:
        qscore = ord(nn_qual) - 33
        
        if qscore > qual_cutoff:
            ngood += 1
    
    qual_perc = round(float(ngood) / read_len, 3)
    return qual_perc

def extract_fasta(sam_fn, qual_perc_cutoff):
    seq_data = {}
    seqid_data = {}
    candidate_fa_fn = sam_fn.replace(".sam", ".candidate.fa")
    candidate_fa_fh = open(candidate_fa_fn, 'w')
    
    with open(sam_fn, 'r') as sam_fh:
        for line in sam_fh:
            lst = line.strip().split('\t')
            seqid = lst[0] + "_" + lst[1]
            chro = lst[2]
            pos = lst[3]
            map_qual = lst[4]
            cigar = lst[5]
            paired_read_chr = lst[6]
            paired_read_pos = lst[7]
            segment_length = lst[8]
            seq = lst[9].upper()
            qual_str = lst[10]
            ngood = 0
            
            if re.findall("^\d+M$", cigar):
                continue
            
            Q20_qual_perc = get_read_qual(qual_str, 20)
            Q30_qual_perc = get_read_qual(qual_str, 30)
            
            if Q20_qual_perc < qual_perc_cutoff:
                continue
            
            if seq not in seq_data:
                seq_data[seq] = {}
                candidate_fa_fh.write(">" + seqid + "\n")
                candidate_fa_fh.write(seq + "\n")
                seqid_data[seqid] = seq
            
            seq_data[seq][seqid] = {"20": Q20_qual_perc, "30": Q30_qual_perc}
    
    candidate_fa_fh.close()
    return seq_data, seqid_data, candidate_fa_fn

def do_blast_indel(WT_ref_fn, mut_ref_fn, candidate_fa_fn):
    WT_blast_fn = WT_ref_fn.replace(".fa", ".blast")
    mut_blast_fn = mut_ref_fn.replace(".fa", ".blast")
    
    command_WT = f"{blast_exe_folder}blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 -gapopen 6 -gapextend 2 -dust no -soft_masking false -query {candidate_fa_fn} -subject {WT_ref_fn} -outfmt 15 > {WT_blast_fn}"
    command_mut = f"{blast_exe_folder}blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 -gapopen 6 -gapextend 2 -dust no -soft_masking false -query {candidate_fa_fn} -subject {mut_ref_fn} -outfmt 15 > {mut_blast_fn}"
    
    os.system(command_WT)
    os.system(command_mut)
    
    return WT_blast_fn, mut_blast_fn

def adjust_indel_pos(qseq, hseq, midline):
    hseq_lst = list(hseq)
    qseq_lst = list(qseq)
    midline_lst = list(midline)
    
    for i in range(len(qseq_lst)):
        if qseq_lst[i] == hseq_lst[i]:
            continue
        
        if hseq_lst[i] == "-" and i < len(hseq_lst) - 1:
            j = 1
            
            while hseq_lst[i + j] == "-" and (i + j) < len(hseq_lst) - 1:
                j += 1
            
            if (i + j) <= (len(hseq_lst) - 1) and hseq_lst[i + j] == qseq_lst[i]:
                hseq_lst[i] = hseq_lst[i + j]
                hseq_lst[i + j] = "-"
                midline_lst[i] = "|"
                midline_lst[i + j] = " "
            else:
                break
        
        if qseq_lst[i] == "-" and (i < len(qseq_lst) - 1):
            j = 1
            
            while (qseq_lst[i + j] == "-") and (i + j) < len(qseq_lst) - 1:
                j += 1
            
            if (i + j) <= (len(qseq_lst) - 1) and qseq_lst[i + j] == hseq_lst[i]:
                qseq_lst[i] = qseq_lst[i + j]
                qseq_lst[i + j] = "-"
                midline_lst[i] = "|"
                midline_lst[i + j] = " "
            else:
                break
    
    adjust_qseq = "".join(qseq_lst)
    adjust_hseq = "".join(hseq_lst)
    adjust_midline = "".join(midline_lst)
    
    return adjust_qseq, adjust_hseq, adjust_midline

def adjust_1seg_indel_pos(seqid, qseq, hseq, midline, geno_pos_inalign, indel, blast_trim = 8):
    hseq_lst = list(hseq)
    qseq_lst = list(qseq)
    midline_lst = list(midline)
    
    for i in range(len(qseq_lst)):
        if qseq_lst[i] == hseq_lst[i]:
            continue
        
        if hseq_lst[i] == "-" and i < len(hseq_lst) - 1:
            j = 1
            
            while hseq_lst[i + j] == "-" and (i + j) < len(hseq_lst) - 1:
                j += 1
            
            if (i + j) <= (len(hseq_lst) - 1) and hseq_lst[i + j] == qseq_lst[i]:
                hseq_lst[i] = hseq_lst[i + j]
                hseq_lst[i + j] = "-"
                midline_lst[i] = "|"
                midline_lst[i + j] = " "
            else:
                break
        
        if qseq_lst[i] == "-" and (i < len(qseq_lst) - 1):
            j = 1
            
            while (qseq_lst[i + j] == "-") and (i + j) < len(qseq_lst) - 1:
                j += 1
            
            if (i + j) <= (len(qseq_lst) - 1) and qseq_lst[i + j] == hseq_lst[i]:
                qseq_lst[i] = qseq_lst[i + j]
                qseq_lst[i + j] = "-"
                midline_lst[i] = "|"
                midline_lst[i + j] = " "
            else:
                break
    
    qseq = "".join(qseq_lst)
    hseq = "".join(hseq_lst)
    midline = "".join(midline_lst)
    
    left_index = 0
    
    for ii in range(len(qseq_lst)):
        if qseq_lst[ii] == hseq_lst[ii]:
            continue
        
        if not (hseq_lst[ii] == "-" or qseq_lst[ii] == "-"):
            continue
        
        left_index = ii - 1
        break
    
    left_qseq = "".join(qseq_lst[0:(left_index + 1)])
    left_hseq = "".join(hseq_lst[0:(left_index + 1)])
    left_midline = "".join(midline_lst[0:(left_index + 1)])
    
    while (left_index - blast_trim + 1 >= 0) and (left_midline[(left_index - blast_trim + 1):(left_index + 1)] != "|" * blast_trim):
        left_index = left_index - 1
    
    indel_left_genome_pos_in_align = []
    left_breakpoint_geno_pos = ".".join(indel.split(".")[0:2])
    j = 0
    
    for bb in hseq:
        if bb == "-":
            indel_left_genome_pos_in_align.append("-")
            continue
        indel_left_genome_pos_in_align.append(geno_pos_inalign[j])
        j += 1
    
    try:
        left_breakpoint_index_inalign = indel_left_genome_pos_in_align.index(left_breakpoint_geno_pos)
        
        if left_breakpoint_index_inalign < left_index:
            left_index = left_breakpoint_index_inalign
        
        left_qseq = "".join(qseq_lst[0:(left_index + 1)])
        left_hseq = "".join(hseq_lst[0:(left_index + 1)])
        left_midline = "".join(midline_lst[0:(left_index + 1)])
        
        right_index = left_index + 1
        right_qseq_rev, right_hseq_rev, right_midline_rev = adjust_indel_pos(qseq[right_index:][::-1], hseq[right_index:][::-1], midline[right_index:][::-1])
        right_qseq, right_hseq, right_midline = right_qseq_rev[::-1], right_hseq_rev[::-1], right_midline_rev[::-1]
        
        adjust_qseq = left_qseq + right_qseq
        adjust_hseq = left_hseq + right_hseq
        adjust_midline = left_midline + right_midline
        keyseq_left_index = left_breakpoint_index_inalign - 6
        keyseq_right_index = left_breakpoint_index_inalign + len(indel.split(".")[-1]) + 6
        key_qseq = adjust_qseq[keyseq_left_index:keyseq_right_index]
        key_hseq = adjust_hseq[keyseq_left_index:keyseq_right_index]
    except:
        adjust_qseq, adjust_hseq, adjust_midline, key_qseq, key_hseq = ["ND", "ND", "ND", "ND", "ND"]
    
    return adjust_qseq, adjust_hseq, adjust_midline, key_qseq, key_hseq

def parse_blast_sv(blast_fn, SV):
    strand_transform={"Plus":"+","Minus":"-","ND":"ND"}
    blast_json=json.load(open(blast_fn))["BlastOutput2"]
    blast_data={}
    for i in range(0,len(blast_json)):
        blast_report= blast_json[i]["report"]
        blast_result= blast_report["results"]
        bl2seq=blast_result["bl2seq"][0]
        seqid=bl2seq["query_title"]
        if seqid not in blast_data:
            blast_data[seqid]={}
        query_len=bl2seq["query_len"]

        hsps="ND"
        align_len="0"
        identity="0"
        gaps="ND"
        query_from="ND"
        query_to="ND"
        query_strand="ND"
        hit_from="ND"
        hit_to="ND"
        hit_strand="ND"
        qseq="ND"
        hseq="ND"
        midline="ND"
        bit_score="ND"
        score="0"
        evalue="ND"

        if ("hits" in bl2seq) and  len(bl2seq["hits"])>0:
            hits=bl2seq["hits"]
            hsps=hits[0]["hsps"][0]
            align_len=hsps["align_len"]
            identity=hsps["identity"]
            gaps=hsps["gaps"]
            query_from=hsps["query_from"]
            query_to=hsps["query_to"]
            query_strand=hsps["query_strand"]
            hit_from=hsps["hit_from"]
            hit_to=hsps["hit_to"]
            hit_strand=hsps["hit_strand"]
            qseq=hsps["qseq"]
            hseq=hsps["hseq"]
            midline=hsps["midline"]
            bit_score=hsps["bit_score"]
            score=hsps["score"]
            evalue=hsps["evalue"]
            
        blast_data[seqid]={
            "initial_SV":SV,
            "align_len":align_len,
            "identity":identity,
            "gaps":gaps,
            "query_from":query_from,
            "query_to":query_to,
            "query_strand":strand_transform[query_strand],
            "hit_from":hit_from,
            "hit_to":hit_to,
            "hit_strand":strand_transform[hit_strand],
            "qseq":qseq,
            "hseq":hseq,
            "midline":midline,
            "bit_score":bit_score,
            "score":score,
            "evalue":evalue
        }
    return blast_data

def parse_blast_indel(blast_fn, indel):
    strand_transform = {"Plus": "+", "Minus": "-", "ND": "ND"}
    blast_json = json.load(open(blast_fn))["BlastOutput2"]
    blast_data = {}
    
    for i in range(0, len(blast_json)):
        blast_report = blast_json[i]["report"]
        blast_result = blast_report["results"]
        bl2seq = blast_result["bl2seq"][0]
        seqid = bl2seq["query_title"]
        
        if seqid not in blast_data:
            blast_data[seqid] = {}
        
        query_len = bl2seq["query_len"]
        hsps = bl2seq.get("hits", [{"hsps": [{}]}])[0]["hsps"][0]
        
        align_len = hsps.get("align_len", "0")
        identity = hsps.get("identity", "0")
        gaps = hsps.get("gaps", "ND")
        query_from = hsps.get("query_from", "ND")
        query_to = hsps.get("query_to", "ND")
        query_strand = strand_transform.get(hsps.get("query_strand"), "ND")
        hit_from = hsps.get("hit_from", "ND")
        hit_to = hsps.get("hit_to", "ND")
        hit_strand = strand_transform.get(hsps.get("hit_strand"), "ND")
        qseq = hsps.get("qseq", "ND")
        hseq = hsps.get("hseq", "ND")
        midline = hsps.get("midline", "ND")
        bit_score = hsps.get("bit_score", "ND")
        score = hsps.get("score", "0")
        evalue = hsps.get("evalue", "ND")
        
        blast_data[seqid] = {
            "initial_indel": indel,
            "align_len": align_len,
            "identity": identity,
            "gaps": gaps,
            "query_from": query_from,
            "query_to": query_to,
            "query_strand": query_strand,
            "hit_from": hit_from,
            "hit_to": hit_to,
            "hit_strand": hit_strand,
            "qseq": qseq,
            "hseq": hseq,
            "midline": midline,
            "bit_score": bit_score,
            "score": score,
            "evalue": evalue
        }
    
    return blast_data

def genotyping_indel(WT_blast_data, mut_blast_data, seqid_data, seq_data, ref_data, indel, output_qual_pairs):
    WT_ref_fn = ref_data["WT_ref_fn"]
    mut_ref_fn = ref_data["mut_ref_fn"]
    WT_ref_seq = ref_data["WT_ref_seq"]
    mut_ref_seq = ref_data["mut_ref_seq"]
    WT_genome_pos = ref_data["WT_genome_pos"]
    mut_genome_pos = ref_data["mut_genome_pos"]
    WT_keyseq = ref_data["WT_keyseq"]
    mut_keyseq = ref_data["mut_keyseq"]
    WT_keyseq_pos = ref_data["WT_keyseq_pos"]
    mut_keyseq_pos = ref_data["mut_keyseq_pos"]
    revcmp_WT_ref_seq = get_revcmp(WT_ref_seq)
    revcmp_mut_ref_seq = get_revcmp(mut_ref_seq)
    revcmp_WT_keyseq = get_revcmp(WT_keyseq)
    revcmp_mut_keyseq = get_revcmp(mut_keyseq)
    indel_read_ids = {"exact": [], "fuzzy": []}
    indel_read_num = {}
    
    for output_qual_pair in output_qual_pairs:
        indel_read_num["exact_" + output_qual_pair] = {}
        indel_read_num["fuzzy_" + output_qual_pair] = {}
    
    for seqid in seqid_data:
        is_indel = "no"
        seq = seqid_data[seqid]
        all_seqids = seq_data[seq]
        WT_blast_info = WT_blast_data[seqid]
        mut_blast_info = mut_blast_data[seqid]
        WT_score = WT_blast_info["score"]
        mut_score = mut_blast_info["score"]
        WT_align_len = WT_blast_info["align_len"]
        mut_align_len = mut_blast_info["align_len"]
        WT_gaps = WT_blast_info["gaps"]
        mut_gaps = mut_blast_info["gaps"]
        WT_qseq = WT_blast_info["qseq"]
        mut_qseq = mut_blast_info["qseq"]
        WT_midline = WT_blast_info["midline"]
        mut_midline = mut_blast_info["midline"]
        WT_hseq = WT_blast_info["hseq"]
        mut_hseq = mut_blast_info["hseq"]
        WT_identity = mut_blast_info["identity"]
        mut_identity = mut_blast_info["identity"]
        WT_hit_strand = WT_blast_info["hit_strand"]
        mut_hit_strand = mut_blast_info["hit_strand"]
        mut_query_from = mut_blast_info["query_from"]
        mut_query_to = mut_blast_info["query_to"]
        mut_hit_from = mut_blast_info["hit_from"]
        mut_hit_to = mut_blast_info["hit_to"]
        WT_query_from = WT_blast_info["query_from"]
        WT_query_to = WT_blast_info["query_to"]

        seqid_score_lst = [(int(WT_score), int(WT_identity), "WT"), (int(mut_score), int(mut_identity), "indel")]
        sorted_seqid_score_lst = sorted(seqid_score_lst, key=lambda x: x[0], reverse=True)
        
        if sorted_seqid_score_lst[0][2] != "indel":
            continue
        
        if sorted_seqid_score_lst[0][1] < len(seq) * 0.95:
            continue
        
        if mut_hit_strand == "+":
            mut_genomepos_inalign = mut_genome_pos[(mut_hit_from - 1):mut_hit_to]
            
            if not (mut_keyseq in mut_qseq.replace("-", "") and mut_keyseq in mut_hseq.replace("-", "")):
                continue
            
            mut_qseq, mut_hseq, mut_midline, key_qseq_in_align, key_hseq_in_align = adjust_1seg_indel_pos(seqid, mut_qseq, mut_hseq, mut_midline, mut_genomepos_inalign, indel)
            
            if "ND" in mut_qseq:
                continue
            
            if mut_keyseq == key_hseq_in_align == key_qseq_in_align:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[0:-1])
                            indel_read_num["exact_" + output_qual_pair][unique_seqidi] = 1
                    
                    indel_read_ids["exact"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
            else:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[0:-1])
                            indel_read_num["fuzzy_" + output_qual_pair][unique_seqidi] = 1
                    
                    indel_read_ids["fuzzy"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
        if mut_hit_strand == "-":
            updated_mut_qseq = get_revcmp(mut_qseq)
            updated_mut_midline = mut_midline[::-1]
            updated_mut_hseq = get_revcmp(mut_hseq)
            updated_mut_query_from = mut_query_to
            updated_mut_query_to = mut_query_from
            updated_mut_hit_from = mut_hit_to
            updated_mut_hit_to = mut_hit_from
            mut_genomepos_inalign = mut_genome_pos[(updated_mut_hit_from - 1):updated_mut_hit_to]
            
            if not (mut_keyseq in updated_mut_qseq.replace("-", "") and mut_keyseq in updated_mut_hseq.replace("-", "")):
                continue
            
            mut_qseq, updated_mut_hseq, updated_mut_midline, key_qseq_in_align, key_hseq_in_align = adjust_1seg_indel_pos(seqid, updated_mut_qseq, updated_mut_hseq, updated_mut_midline, mut_genomepos_inalign, indel)
            
            if "ND" in mut_qseq:
                continue
            
            if mut_keyseq == key_hseq_in_align == key_qseq_in_align:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[0:-1])
                            indel_read_num["exact_" + output_qual_pair][unique_seqidi] = 1
                    
                    indel_read_ids["exact"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
            else:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[0:-1])
                            indel_read_num["fuzzy_" + output_qual_pair][unique_seqidi] = 1
                    
                    indel_read_ids["fuzzy"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
    
    exact_indel_read_ids = indel_read_ids["exact"]
    fuzzy_indel_read_ids = indel_read_ids["fuzzy"]
    exact_indel_read_num = [len(indel_read_num["exact_" + op]) for op in output_qual_pairs]
    fuzzy_indel_read_num = [(len(indel_read_num["fuzzy_" + op]) + len(indel_read_num["exact_" + op])) for op in output_qual_pairs]
    
    return (exact_indel_read_ids, fuzzy_indel_read_ids, exact_indel_read_num, fuzzy_indel_read_num)


def process_indel(indel_line, output_qual_pairs, centromere_regions, clean_deepseq_folder, qual_perc_cutoff, out_indel_reads_fh):
    lst = indel_line.strip().split("\t")
    initial_indel = lst[0]
    indel = lst[1]
    indel_padding = lst[2]
    bam = lst[3]
    sample = os.path.basename(bam).replace(".bam", "")
    bit_fn, genome_version, chr_in_bam = get_bit_fn(bam)
    chro, pos, ref_allele, mut_allele = indel.split(".")

    out_indel = (f"{sample}\t{initial_indel}\t{indel}\t{indel_padding}\tND"+ "\tND" * len(output_qual_pairs)+ "\tND\n")

    if chro not in centromere_regions[genome_version] or (
        centromere_regions[genome_version][chro][0] < int(pos) < centromere_regions[genome_version][chro][1]
    ):
        return out_indel

    sam_fn = extract_sam_indel(indel, bam, chr_in_bam)
    depth_pos = get_depth_indel(indel, sample, chr_in_bam, bam, clean_deepseq_folder)
    seq_data, seqid_data, candidate_fa_fn = extract_fasta(sam_fn, qual_perc_cutoff)
    candidate_fa_fn_check_command = f"wc -l {candidate_fa_fn}"
    candidate_fa_fn_check_result = str(
        subprocess.check_output(candidate_fa_fn_check_command, shell=True), "utf-8"
    )
    files_to_delete = [sam_fn, candidate_fa_fn]

    if candidate_fa_fn_check_result.split(" ")[0] == "0":
        return (f"{out_indel}{sample}\t{initial_indel}\t{indel}\t{indel_padding}\t0"+ "\t0" * len(output_qual_pairs)+ f"\t{depth_pos}\n")

    print("   #reads: ", candidate_fa_fn_check_result.split(" ")[0])

    ref_data = make_ref_indel(indel, bam, bit_fn)
    WT_ref_fn = ref_data["WT_ref_fn"]
    mut_ref_fn = ref_data["mut_ref_fn"]

    WT_blast_fn, mut_blast_fn = do_blast_indel(WT_ref_fn, mut_ref_fn, candidate_fa_fn)
    files_to_delete.extend([WT_ref_fn, mut_ref_fn, WT_blast_fn, mut_blast_fn])

    WT_blast_data = parse_blast_indel(WT_blast_fn, indel)
    mut_blast_data = parse_blast_indel(mut_blast_fn, indel)
    (
        exact_indel_read_ids,
        fuzzy_indel_read_ids,
        exact_indel_read_num,
        fuzzy_indel_read_num,
    ) = genotyping_indel(
        WT_blast_data, mut_blast_data, seqid_data, seq_data, ref_data, indel, output_qual_pairs
    )

    exact_indel_read_num_str = '\t'.join(map(str, exact_indel_read_num))
    fuzzy_indel_read_num_str = '\t'.join(map(str, fuzzy_indel_read_num))
    out_indel = (
        sample + "\t" + initial_indel + "\t" + indel + "\t" + indel_padding + "\t" +
        exact_indel_read_num_str + "\t" + fuzzy_indel_read_num_str + "\t" + depth_pos + "\n"
    )

    out_indel_reads = ""
    for read_id in exact_indel_read_ids:
        out_indel_reads_fh.write(initial_indel+"\t"+indel+"\t"+indel_padding+"\t"+("\t").join(read_id)+"\t"+bam+"\texact_match\n")
    for read_id in fuzzy_indel_read_ids:
        out_indel_reads_fh.write(initial_indel+"\t"+indel+"\t"+indel_padding+"\t"+("\t").join(read_id)+"\t"+bam+"\tfuzzy_match\n")

    # Remove files
    for file_to_delete in set(files_to_delete):
        if os.path.isfile(file_to_delete):
            os.remove(file_to_delete)

    return f"{out_indel}{out_indel_reads}"

def adjust_SV(SV):
    if SV.count(".") == 5:
        SV += "."
    cmp_strand = {"+": "-", "-": "+"}
    SV_lst = SV.split(".")
    chrA = SV_lst[0]
    if "chr" not in chrA:
        chrA = "chr" + chrA
    if chrA == "chrX":
        chrA = "chr23"
    if chrA == "chrY":
        chrA = "chr24"
    posA = SV_lst[1]
    strandA = SV_lst[2]
    chrB = SV_lst[3]
    if "chr" not in chrB:
        chrB = "chr" + chrB
    if chrB == "chrX":
        chrB = "chr23"
    if chrB == "chrY":
        chrB = "chr24"
    posB = SV_lst[4]
    strandB = SV_lst[5]
    insertion = SV_lst[6]

    adjusted_chrA = chrA
    adjusted_posA = posA
    adjusted_strandA = strandA
    adjusted_chrB = chrB
    adjusted_posB = posB
    adjusted_strandB = strandB
    adjusted_insertion = insertion

    if (
        int(chrA.replace("chr", "")) > int(chrB.replace("chr", ""))
        or (
            int(chrA.replace("chr", ""))
            == int(chrB.replace("chr", ""))
            and int(posA) > int(posB)
        )
    ):
        adjusted_chrA = chrB
        adjusted_posA = posB
        adjusted_strandA = cmp_strand[strandB]
        adjusted_chrB = chrA
        adjusted_posB = posA
        adjusted_strandB = cmp_strand[strandA]
        adjusted_insertion = get_revcmp(insertion)
    if adjusted_chrA == "chr23":
        adjusted_chrA = "chrX"
    if adjusted_chrA == "chr24":
        adjusted_chrA = "chrY"
    if adjusted_chrB == "chr23":
        adjusted_chrB = "chrX"
    if adjusted_chrB == "chr24":
        adjusted_chrB = "chrY"

    adjusted_SV = f"{adjusted_chrA}.{adjusted_posA}.{adjusted_strandA}.{adjusted_chrB}.{adjusted_posB}.{adjusted_strandB}.{adjusted_insertion}"
    return adjusted_SV

def get_depth_sv(SV, sample, chr_in_bam, bam, clean_deepseq_folder="no"):
    chrA, posA, strandA, chrB, posB, strandB, insertion = SV.split(".")
    depth_posA = "0"
    depth_posB = "0"
    
    def get_depth_command(chr, pos, bam):
        return (
            f"samtools depth -r {chr.replace('chr', '')}:{pos}-{pos} {bam}|cut -f3"
        )

    def get_tabix_command(file_path, chr, pos):
        return f"tabix {file_path} {chr.replace('chr', '')}:{pos}-{pos}"

    if not clean_deepseq_folder == "no":
        clean_deepseq_fn = f"{clean_deepseq_folder}/{sample}.count.gz"
        
        try:
            posA_result = (
                str(subprocess.check_output(get_tabix_command(clean_deepseq_fn, chrA, posA), shell=True), "utf-8").strip()
            )
            if len(posA_result) > 0:
                depth_posA = posA_result.split("\n")[0].split("\t")[-1]
        except:
            print(f"no depth: {chrA}.{posA}")

        try:
            posB_result = (
                str(subprocess.check_output(get_tabix_command(clean_deepseq_fn, chrB, posB), shell=True), "utf-8").strip()
            )
            if len(posB_result) > 0:
                depth_posB = posB_result.split("\n")[0].split("\t")[-1]
        except:
            print(f"no depth: {chrB}.{posB}")
    else:
        depth_posA_command = get_depth_command(chrA, posA, bam)
        depth_posB_command = get_depth_command(chrB, posB, bam)

        try:
            depth_posA = str(subprocess.check_output(depth_posA_command, shell=True), "utf-8").strip()
        except:
            print(f"no depth: {chrA}.{posA}")

        try:
            depth_posB = str(subprocess.check_output(depth_posB_command, shell=True), "utf-8").strip()
        except:
            print(f"no depth: {chrB}.{posB}")
    
    if depth_posA == "":
        depth_posA = "0"
    if depth_posB == "":
        depth_posB = "0"
        
    return depth_posA, depth_posB

def make_ref_sv(SV, bam, bit_fn):
    SV_lst = SV.split(".")
    chrA, posA, strandA, chrB, posB, strandB, insertion = SV_lst
    sample = os.path.basename(bam)
    if "chr" not in SV:
        chrA = "chr" + chrA
        chrB = "chr" + chrB
    startA = int(posA) - 200
    endA = int(posA) + 200
    startB = int(posB) - 200
    endB = int(posB) + 200
    posA_ref_command = f"samtools faidx {bit_fn} {chrA}:{startA}-{endA}"
    posB_ref_command = f"samtools faidx {bit_fn} {chrB}:{startB}-{endB}"
    posA_ref_seq = ("" .join(str(subprocess.check_output(posA_ref_command, shell=True), "utf-8").strip().split("\n")[1:])).upper()
    posB_ref_seq = ("" .join(str(subprocess.check_output(posB_ref_command, shell=True), "utf-8").strip().split("\n")[1:])).upper()
    posA_genome_pos = [f"{chrA}.{i}" for i in range(startA, endA + 1)]
    posB_genome_pos = [f"{chrB}.{i}" for i in range(startB, endB + 1)]
    insertion_pos = [f"ins.{i}" for i in range(len(insertion))]
    
    if strandA == "+" and strandB == "-":
        posB_ref_seq = get_revcmp(posB_ref_seq)
        posB_genome_pos = posB_genome_pos[::-1]
    if strandA == "-" and strandB == "+":
        posA_ref_seq = get_revcmp(posA_ref_seq)
        posA_genome_pos = posA_genome_pos[::-1]
    if strandA == "-" and strandB == "-":
        posA_ref_seq = get_revcmp(posA_ref_seq)
        posA_genome_pos = posA_genome_pos[::-1]
        posB_ref_seq = get_revcmp(posB_ref_seq)
        posB_genome_pos = posB_genome_pos[::-1]
    
    SV_ref_seq = posA_ref_seq[0:201] + insertion + posB_ref_seq[200:]
    SV_genome_pos = posA_genome_pos[0:201] + insertion_pos + posB_genome_pos[200:]
    
    posA_ref_fn = "".join(random.choice(string.ascii_lowercase) for x in range(7)) + "." + sample + "." + chrA + "." + posA + "." + strandA + "_ref.fa"
    posB_ref_fn = "".join(random.choice(string.ascii_lowercase) for x in range(6)) + "." + sample + "." + chrB + "." + posB + "." + strandB + "_ref.fa"
    SV_ref_fn = "".join(random.choice(string.ascii_lowercase) for x in range(8)) + "." + sample + "." + SV + "_ref.fa"
    
    with open(posA_ref_fn, "w") as posA_ref_fh:
        posA_ref_fh.write(f">{sample}.{chrA}.{posA}.{strandA}\n")
        posA_ref_fh.write(posA_ref_seq + "\n")
    with open(posB_ref_fn, "w") as posB_ref_fh:
        posB_ref_fh.write(f">{sample}.{chrB}.{posB}.{strandB}\n")
        posB_ref_fh.write(posB_ref_seq + "\n")
    with open(SV_ref_fn, "w") as SV_ref_fh:
        SV_ref_fh.write(f">{sample}.{SV}\n")
        SV_ref_fh.write(SV_ref_seq + "\n")
    
    ref_data = {
        "posA_ref_fn": posA_ref_fn,
        "posB_ref_fn": posB_ref_fn,
        "SV_ref_fn": SV_ref_fn,
        "posA_ref_seq": posA_ref_seq,
        "posB_ref_seq": posB_ref_seq,
        "SV_ref_seq": SV_ref_seq,
        "insertion_seq": insertion,
        "posA_genome_pos": posA_genome_pos,
        "posB_genome_pos": posB_genome_pos,
        "SV_genome_pos": SV_genome_pos,
        "insertion_pos": insertion_pos,
        "posA_keyseq": posA_ref_seq[(200 - 6):(200 + 6)],
        "posA_keyseq_pos": posA_genome_pos[(200 - 6):(200 + 6)],
        "posB_keyseq": posB_ref_seq[(200 - 6):(200 + 6)],
        "posB_keyseq_pos": posB_genome_pos[(200 - 6):(200 + 6)],
        "SV_keyseq": SV_ref_seq[(200 - 6):(200 + len(insertion) + 6)],
        "SV_keyseq_pos": SV_genome_pos[(200 - 6):(200 + len(insertion) + 6)]
    }
    
    return ref_data


def extract_sam_sv(SV, bam, chr_in_bam):
    SV_lst = SV.split(".")
    chrA, posA, strandA, chrB, posB, strandB, _ = SV_lst
    sample = os.path.basename(bam).replace(".bam", "")
    
    def get_chromosome(chr_in_bam, chr):
        if chr_in_bam == "yes" and ("chr" not in chr):
            return "chr" + chr
        elif chr_in_bam == "no" and ("chr" in chr):
            return chr.replace("chr", "")
        return chr

    chrA = get_chromosome(chr_in_bam, chrA)
    chrB = get_chromosome(chr_in_bam, chrB)
    
    sam_fn = ('').join(random.choice(string.ascii_lowercase) for x in range(8)) + f".{SV}.{sample}.sam"
    command_A = f"samtools view {bam} {chrA}:{int(posA) - 100}-{int(posA) + 100} | sort -k1,1 > {sam_fn}"
    command_B = f"samtools view {bam} {chrB}:{int(posB) - 100}-{int(posB) + 100} | sort -k1,1 >> {sam_fn}"
    
    os.system(command_A)
    os.system(command_B)
    
    return sam_fn


def do_blast_sv(posA_ref_fn, posB_ref_fn, SV_ref_fn, candidate_fa_fn):
    posA_blast_fn = posA_ref_fn.replace(".fa", ".blast")
    posB_blast_fn = posB_ref_fn.replace(".fa", ".blast")
    SV_blast_fn = SV_ref_fn.replace(".fa", ".blast")
    
    blast_command = (
        "{}blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 "
        "-gapopen 6 -gapextend 2 -dust no -soft_masking false -query {} "
        "-subject {} -outfmt 15 > {}"
    )
    
    os.system(blast_command.format(blast_exe_folder, candidate_fa_fn, posA_ref_fn, posA_blast_fn))
    os.system(blast_command.format(blast_exe_folder, candidate_fa_fn, posB_ref_fn, posB_blast_fn))
    os.system(blast_command.format(blast_exe_folder, candidate_fa_fn, SV_ref_fn, SV_blast_fn))
    
    return (posA_blast_fn, posB_blast_fn, SV_blast_fn)

def genotyping_SV(posA_blast_data, posB_blast_data, SV_blast_data, seqid_data, seq_data, ref_data, SV, output_qual_pairs):
    posA_ref_fn, posB_ref_fn, SV_ref_fn = ref_data["posA_ref_fn"], ref_data["posB_ref_fn"], ref_data["SV_ref_fn"]
    posA_ref_seq, posB_ref_seq, SV_ref_seq = ref_data["posA_ref_seq"], ref_data["posB_ref_seq"], ref_data["SV_ref_seq"]
    insertion_seq, posA_genome_pos, posB_genome_pos = ref_data["insertion_seq"], ref_data["posA_genome_pos"], ref_data["posB_genome_pos"]
    SV_genome_pos, insertion_pos = ref_data["SV_genome_pos"], ref_data["posA_genome_pos"]
    posA_keyseq, posB_keyseq, SV_keyseq = ref_data["posA_keyseq"], ref_data["posB_keyseq"], ref_data["SV_keyseq"]
    posA_keyseq_pos, posB_keyseq_pos, SV_keyseq_pos = ref_data["posA_keyseq_pos"], ref_data["posB_keyseq_pos"], ref_data["SV_keyseq_pos"]

    SV_read_ids = {"exact": [], "fuzzy": []}
    SV_read_num = {}

    for output_qual_pair in output_qual_pairs:
        SV_read_num[f"exact_{output_qual_pair}"] = {}
        SV_read_num[f"fuzzy_{output_qual_pair}"] = {}

    for seqid in seqid_data:
        seq = seqid_data[seqid]
        all_seqids = seq_data[seq]

        posA_blast_info, posB_blast_info, SV_blast_info = posA_blast_data[seqid], posB_blast_data[seqid], SV_blast_data[seqid]
        posA_score, posB_score, SV_score = posA_blast_info["score"], posB_blast_info["score"], SV_blast_info["score"]
        posA_identity, posB_identity, SV_identity = posA_blast_info["identity"], posB_blast_info["identity"], SV_blast_info["identity"]
        posA_qseq, posB_qseq, SV_qseq = posA_blast_info["qseq"], posB_blast_info["qseq"], SV_blast_info["qseq"]
        posA_midline, posB_midline, SV_midline = posA_blast_info["midline"], posB_blast_info["midline"], SV_blast_info["midline"]
        posA_hseq, posB_hseq, SV_hseq = posA_blast_info["hseq"], posB_blast_info["hseq"], SV_blast_info["hseq"]
        posA_hit_strand, posB_hit_strand, SV_hit_strand = posA_blast_info["hit_strand"], posB_blast_info["hit_strand"], SV_blast_info["hit_strand"]
        
        SV_query_from, SV_query_to, SV_hit_from, SV_hit_to = SV_blast_info["query_from"], SV_blast_info["query_to"], SV_blast_info["hit_from"], SV_blast_info["hit_to"]

        seqid_score_lst = [(int(posA_score), int(posA_identity), "posA"), (int(posB_score), int(posB_identity), "posB"), (int(SV_score), int(SV_identity), "SV")]
        sorted_seqid_score_lst = sorted(seqid_score_lst, key=lambda x: x[0], reverse=True)

        if sorted_seqid_score_lst[0][2] != "SV":
            continue

        if sorted_seqid_score_lst[0][1] < len(seq) * 0.95:
            continue

        SV_genomepos_inalign = SV_genome_pos[(SV_hit_from - 1):SV_hit_to]

        if SV_hit_strand == "+":
            if SV_keyseq not in SV_qseq.replace("-", "") or SV_keyseq not in SV_hseq.replace("-", ""):
                continue
            SV_qseq, SV_hseq, SV_midline, key_qseq_in_align, key_hseq_in_align = adjust_1seg_indel_pos(seqid, SV_qseq, SV_hseq, SV_midline, SV_genomepos_inalign, SV)
            if "ND" in SV_qseq:
                continue
            if SV_keyseq == key_hseq_in_align == key_qseq_in_align:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[:-1])
                            SV_read_num[f"exact_{output_qual_pair}"][unique_seqidi] = 1
                    SV_read_ids["exact"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
            else:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[:-1])
                            SV_read_num[f"fuzzy_{output_qual_pair}"][unique_seqidi] = 1
                    SV_read_ids["fuzzy"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
        elif SV_hit_strand == "-":
            updated_SV_qseq = get_revcmp(SV_qseq)
            updated_SV_midline = SV_midline[::-1]
            updated_SV_hseq = get_revcmp(SV_hseq)
            updated_SV_query_from = SV_query_to
            updated_SV_query_to = SV_query_from
            updated_SV_hit_from = SV_hit_to
            updated_SV_hit_to = SV_hit_from
            SV_genomepos_inalign = SV_genome_pos[(updated_SV_hit_from - 1):updated_SV_hit_to]

            if SV_keyseq not in updated_SV_qseq.replace("-", "") or SV_keyseq not in updated_SV_hseq.replace("-", ""):
                continue

            SV_qseq, updated_SV_hseq, updated_SV_midline, key_qseq_in_align, key_hseq_in_align = adjust_1seg_indel_pos(seqid, updated_SV_qseq, updated_SV_hseq, updated_SV_midline, SV_genomepos_inalign, SV)
            if "ND" in SV_qseq:
                continue

            if SV_keyseq == key_hseq_in_align == key_qseq_in_align:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[:-1])
                            SV_read_num[f"exact_{output_qual_pair}"][unique_seqidi] = 1
                    SV_read_ids["exact"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])
            else:
                for seqidi in all_seqids:
                    for output_qual_pair in output_qual_pairs:
                        seqidi_qual_cutoff, seqidi_perc_cutoff = output_qual_pair.split("_")
                        if all_seqids[seqidi][seqidi_qual_cutoff] > float(seqidi_perc_cutoff):
                            unique_seqidi = "_".join(seqidi.split("_")[:-1])
                            SV_read_num[f"fuzzy_{output_qual_pair}"][unique_seqidi] = 1
                    SV_read_ids["fuzzy"].append([seqidi, str(all_seqids[seqidi]["20"]), str(all_seqids[seqidi]["30"])])

    exact_SV_read_ids = SV_read_ids["exact"]
    fuzzy_SV_read_ids = SV_read_ids["fuzzy"]
    exact_SV_read_num = [len(SV_read_num[f"exact_{op}"]) for op in output_qual_pairs]
    fuzzy_SV_read_num = [(len(SV_read_num[f"fuzzy_{op}"]) + len(SV_read_num[f"exact_{op}"])) for op in output_qual_pairs]

    return exact_SV_read_ids, fuzzy_SV_read_ids, exact_SV_read_num, fuzzy_SV_read_num

def process_sv(line, centromere_regions, output_qual_pair_SV_num, out_SV_fh, out_SV_reads_fh, qual_perc_cutoff, clean_deepseq_folder, output_qual_pairs):
    files_to_delete = []
    lst = line.strip().split("\t")
    initial_SV = lst[0]
    SV = lst[1]
    bam = lst[2]
    sample = os.path.basename(bam).replace(".bam", "")
    bit_fn, genome_version, chr_in_bam = get_bit_fn(bam)

    if initial_SV.count(".") == 5:
        initial_SV = initial_SV + "."

    chrA, posA, strandA, chrB, posB, strandB, insertion = SV.split(".")
    
    if not (chrA in centromere_regions[genome_version] and chrB in centromere_regions[genome_version]):
        out_SV = f"{sample}\t{initial_SV}\t{SV}\tND" * len(output_qual_pair_SV_num) + "\tND\tND"
        out_SV_fh.write(out_SV + "\n")
        return

    posA_range = centromere_regions[genome_version][chrA]
    posB_range = centromere_regions[genome_version][chrB]
    if posA_range[0] < int(posA) < posA_range[1] or posB_range[0] < int(posB) < posB_range[1]:
        out_SV = f"{sample}\t{initial_SV}\t{SV}\tND" * len(output_qual_pair_SV_num) + "\tND\tND"
        out_SV_fh.write(out_SV + "\n")
        return
    
    sam_fn = extract_sam_sv(SV, bam, chr_in_bam)
    sam_fn_check_command="wc -l "+sam_fn
    sam_fn_check_result=str(subprocess.check_output(sam_fn_check_command,shell=True),"utf-8")

    print("   #reads: ", sam_fn_check_result.split(" ")[0])

    seq_data, seqid_data, candidate_fa_fn = extract_fasta(sam_fn, qual_perc_cutoff)
    files_to_delete = [sam_fn, candidate_fa_fn]
    candidate_fa_fn_check_command = "wc -l " + candidate_fa_fn
    candidate_fa_fn_check_result = str(subprocess.check_output(candidate_fa_fn_check_command, shell=True), "utf-8")
    depth_posA, depth_posB = get_depth_sv(SV, sample, chr_in_bam, bam, clean_deepseq_folder)

    if candidate_fa_fn_check_result.split(" ")[0] == "0":
        SV_reads = "0"
        out_SV = f"{sample}\t{initial_SV}\t{SV}\t0" * len(output_qual_pair_SV_num) + f"\t{depth_posA}\t{depth_posB}\n"
        out_SV_fh.write(out_SV)
        os.remove(candidate_fa_fn)
        os.remove(sam_fn)
        return

    ref_data = make_ref_sv(SV, bam, bit_fn)
    posA_ref_fn = ref_data["posA_ref_fn"]
    posB_ref_fn = ref_data["posB_ref_fn"]
    SV_ref_fn = ref_data["SV_ref_fn"]

    posA_blast_fn, posB_blast_fn, SV_blast_fn = do_blast_sv(posA_ref_fn, posB_ref_fn, SV_ref_fn, candidate_fa_fn)
    files_to_delete += [posA_ref_fn, posB_ref_fn, SV_ref_fn, posA_blast_fn, posB_blast_fn, SV_blast_fn]

    posA_blast_data = parse_blast_sv(posA_blast_fn, SV)
    posB_blast_data = parse_blast_sv(posB_blast_fn, SV)
    SV_blast_data = parse_blast_sv(SV_blast_fn, SV)
    exact_SV_read_ids, fuzzy_SV_read_ids, exact_SV_read_num, fuzzy_SV_read_num = genotyping_SV(
        posA_blast_data, posB_blast_data, SV_blast_data, seqid_data, seq_data, ref_data, SV, output_qual_pairs
    )

    out_SV = sample + "\t" + initial_SV + "\t" + SV + "\t" + "\t".join(map(str, exact_SV_read_num)) + "\t" + "\t".join(map(str, fuzzy_SV_read_num)) + "\t" + depth_posA + "\t" + depth_posB + "\n"

    out_SV_fh.write(out_SV)
    for read_id in exact_SV_read_ids:
        out_SV_reads_fh.write("{}\t{}\t{}\t{}\texact_match\n".format(initial_SV, SV, '\t'.join(read_id), bam))

    for read_id in fuzzy_SV_read_ids:
        out_SV_reads_fh.write("{}\t{}\t{}\t{}\tfuzzy_match\n".format(initial_SV, SV, '\t'.join(read_id), bam))

    for file_to_delete in set(files_to_delete):
        os.remove(file_to_delete)

def main():
    if (len(sys.argv) > 4) or (len(sys.argv) < 3):
        print("Please input [SV/Indel] [SVIndel list fn] [optional: clean_deepseq_folder]")
        print_usage()

    SVIndel_fn = sys.argv[2]
    clean_deepseq_folder = sys.argv[3] if len(sys.argv) == 4 else "no"

    SVIndel_padding = 60  # Final SV allele with 60 bp flanking seq chrA.posA.starndA.chrB.posB.starndB.flanking_A_seq.insertion_seq.flanking_B_seq
    blast_trim = 8  # Trim blast midlines at breakpoints to 8bp exact match '||||||||'
    flanking_length = 200
    qual_cutoff = 20
    qual_perc_cutoff = 0.8
    output_qual_pairs = ["20_0.8", "20_0.9", "30_0.9", "30_0.95"]

    # Define centromere regions for different reference genomes
    centromere_regions = {
        "hg19": {
            "chr1": [120000000, 145000000],
            "chr2": [87000000, 98000000],
            "chr3": [88000000, 96000000],
            "chr4": [48000000, 52900000],
            "chr5": [44000000, 52000000],
            "chr6": [57000000, 64000000],
            "chr7": [57000000, 64000000],
            "chr8": [42000000, 48000000],
            "chr9": [44000000, 70000000],
            "chr10": [38000000, 44000000],
            "chr11": [50000000, 56000000],
            "chr12": [34000000, 40000000],
            "chr13": [0, 20000000],
            "chr14": [0, 20000000],
            "chr15": [0, 20000000],
            "chr16": [32000000, 49500000],
            "chr17": [22000000, 26000000],
            "chr18": [14000000, 20000000],
            "chr19": [23000000, 29000000],
            "chr20": [25000000, 30000000],
            "chr21": [0, 15000000],
            "chr22": [0, 17000000],
            "chrX": [57000000, 62500000],
            "chrY": [9000000, 14000000]
        },
        "hg38": {
            "chr1": [120000000, 145000000],
            "chr2": [87000000, 98000000],
            "chr3": [88000000, 96000000],
            "chr4": [48600000, 54000000],
            "chr5": [44000000, 52000000],
            "chr6": [57000000, 62000000],
            "chr7": [57000000, 64000000],
            "chr8": [43000000, 46400000],
            "chr9": [42000000, 70000000],
            "chr10": [38000000, 44000000],
            "chr11": [50000000, 56000000],
            "chr12": [32000000, 40000000],
            "chr13": [0, 20000000],
            "chr14": [0, 20000000],
            "chr15": [0, 20000000],
            "chr16": [34000000, 48000000],
            "chr17": [22000000, 28000000],
            "chr18": [14000000, 22000000],
            "chr19": [23000000, 29000000],
            "chr20": [25000000, 33000000],
            "chr21": [0, 14000000],
            "chr22": [0, 17000000],
            "chrX": [57000000, 63000000],
            "chrY": [6000000, 11000000]
        }
    }

    # for Indel geno
    if (sys.argv[1] == 'indel') or (sys.argv[1] == 'Indel') or (sys.argv[1] == 'INDEL'):
        print("+++Indel Genotype for: ", SVIndel_fn)

        out_indel_fn = f"{SVIndel_fn}_genotyping"
        output_qual_pair_indel_num = [
            f"exact_match_{output_qual_pair}" for output_qual_pair in output_qual_pairs
            ] + [
            f"fuzzy_match_{output_qual_pair}" for output_qual_pair in output_qual_pairs
            ]
        out_indel_fh = open(out_indel_fn, "w")
        out_indel_fh.write(
            "sample\tinitial_indel\tcurated_indel\tindel_padding\t"
            + "\t".join(output_qual_pair_indel_num)
            + "\tdepth_pos\n"
        )

        out_indel_reads_fn = f"{SVIndel_fn}_genotyping_indel_reads"
        out_indel_reads_fh = open(out_indel_reads_fn, "w")
        out_indel_reads_fh.write("indel\tindel_padding\tread_id\tQ20_perc\tQ30_perc\tbam\tcomment\n")

        with open(SVIndel_fn, "r") as indel_fh:
            for line in indel_fh:
                processed_output = process_indel(line, output_qual_pairs, centromere_regions, clean_deepseq_folder, qual_perc_cutoff, out_indel_reads_fh)
                if processed_output:
                    out_indel_fh.write(processed_output)

        out_indel_fh.close()
        out_indel_reads_fh.close()

    # for SV geno
    if (sys.argv[1] == 'sv') or (sys.argv[1] == 'Sv') or (sys.argv[1] == 'SV'):
        print("+++SV Genotype for: ", SVIndel_fn)

        out_SV_fn = SVIndel_fn + "_genotyping"
        output_qual_pair_SV_num = ["exact_match_" + output_qual_pair for output_qual_pair in output_qual_pairs] + ["fuzzy_match_" + output_qual_pair for output_qual_pair in output_qual_pairs]
        out_SV_header = "sample\tinitial_SV\tcurated_SV\t" + "\t".join(output_qual_pair_SV_num) + "\tdepth_posA\tdepth_posB"
        out_SV_fh=open(out_SV_fn,"w")
        out_SV_fh.write(out_SV_header+"\n")

        out_SV_reads_fn = SVIndel_fn + "_genotyping_SV_reads"
        out_SV_reads_header = "initial_SV\tSV\tread_id\tQ20_perc\tQ30_perc\tbam\tcomment"
        out_SV_reads_fh=open(out_SV_reads_fn,"w")
        out_SV_reads_fh.write(out_SV_reads_header+"\n")

        with open(SVIndel_fn, "r") as SV_fh:
            for line in SV_fh:
                process_sv(line, centromere_regions, output_qual_pair_SV_num, out_SV_fh, out_SV_reads_fh, qual_perc_cutoff, clean_deepseq_folder, output_qual_pairs)

        out_SV_fh.close()
        out_SV_reads_fh.close()

    print("Genotyping is complete!")

if __name__ == "__main__":
    main()
