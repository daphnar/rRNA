import os
from Bio import SeqIO,Seq,SeqRecord
import glob
import pandas as pd
from datetime import datetime
import subprocess
import sys
def run_shell_command(command,debug=True):
    if debug:
        print('Running: %s'%command)
        print('Starting: %s\n\n'%datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    output = subprocess.check_output(command, shell = True)
    if debug:
        print('End: %s\n\n'%datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    return output

def split_allseqs_to_single_files(all_seq_fasta,output_dir):
    f_in = SeqIO.parse(all_seq_fasta,'fasta')
    output_files = []
    for record in f_in:
        output_file =os.path.join(output_dir,#'tmp',
                record.id.replace('|','___').replace('/','_slash_').replace(';size=','_size_')+'.fa')
        if not os.path.exists(output_file):
            with open(output_file,'w') as f_out:
                SeqIO.write([record],f_out,'fasta')
        output_files.append(output_file)
    return output_files

OAK_BASE=''
base_path=OAK_BASE

ribosome_reference_path=os.path.join(base_path,'shared/genomes/ribosome')
output_dir='output_dir'
msa_path = os.path.join(base_path,MSA')
output_tmp_dir=os.path.join(base_path,output_dir,'tmp'
 
helix_pos_table=os.path.join(base_path,\
    'supp_tables/TableS1-helix-annotations.csv')
es_pos_table=os.path.join(base_path,\
    'supp_tables/TableS2-ES-core-annotations.csv')

def dereplicate_remote():
    all_reference_records={}
    for region_file in ['helix_Theo.fa','es_Theo.fa','non_es_Theo.fa']:
        all_reference_records.update({record.id:record for record in \
            SeqIO.parse(os.path.join(ribosome_reference_path,region_file),'fasta')})
    for region_atlas in glob.glob(os.path.join(atlas_path,'*.fa')):
        region=os.path.basename(region_atlas).replace('.atlas.fa','')
        if region=='all_region':
            continue
        united_variants=[all_reference_records[region]]+list(SeqIO.parse(region_atlas,'fasta'))
        united_output=os.path.join(output_dir,region+'.fa')
        with open(united_output,'w') as f_h:
            SeqIO.write(united_variants,f_h,'fasta')
        dereplicated_output = os.path.join(output_dir, region + '.dereplication.fa')
        command = '/oak/stanford/groups/pritch/users/daphna/tools/vsearch-2.17.1/bin/vsearch --minseqlength 2 --derep_fulllength %s --sizeout --fasta_width 0 --output %s' \
                  % (united_output,dereplicated_output)
        os.system(command)

def remove_imputed_variants():
    for region_atlas in glob.glob(os.path.join(output_dir, '*.dereplication.fa')):
        keep_reads=[]
        for record in SeqIO.parse(region_atlas,'fasta'):
            if record.id.find('imputed')==-1:
                keep_reads.append(record)
        with open(region_atlas, 'w') as f_h:
            SeqIO.write(keep_reads, f_h, 'fasta')

def run_needle_on_region_variants(region_atlas,output_dir,tmp_dir,exe_path=''):
    executable = exe_path + 'needle'
    output_region_alignment = region_atlas.replace('dereplication', 'Needle_SA')
    all_aligned_records = []
    seq_files = split_allseqs_to_single_files(region_atlas, tmp_dir)
    ref_file = seq_files[0]
    other_variants = seq_files[1:]
    all_gap_ref_positions = []
    for idx,variant_file in enumerate(other_variants):
        command = '%s -asequence %s -bsequence %s -datafile EDNAFULL -gapopen 2' \
                  ' -gapextend 10 -snucleotide1' \
                  ' -snucleotide2 -outfile %s -aformat fasta' % (
                  executable, ref_file, variant_file, variant_file + '.needle')
        if not os.path.exists(variant_file + '.needle'):
            run_shell_command(command, False)
        ref_record_s = pd.Series(list(list(SeqIO.parse(variant_file + '.needle', 'fasta'))[0].seq))
        ref_record_gap = list(ref_record_s[ref_record_s == '-'].index.values)
        for gap_n, gap_pos in enumerate(ref_record_gap):
            ref_record_gap[gap_n] = gap_pos - gap_n  # gaps make the sequence longer which we correct
        how_many_times = pd.Series(ref_record_gap).value_counts().astype(int)
        all_gap_ref_positions.append(how_many_times)
        if idx%5000==0:
            print('5000 done')
    try:
        all_gap_ref_positions = pd.concat(all_gap_ref_positions).sort_index()
    except ValueError:
        print(region_atlas)
        return
    max_gap_positions = []
    for pos, pos_df in all_gap_ref_positions.reset_index().groupby('index'):
        max_gap_size_idx = pos_df[0].idxmax()
        max_gap_positions.append([int(pos), int(pos_df.loc[max_gap_size_idx, 0])])
    all_gap_ref_positions_max = pd.DataFrame(max_gap_positions, columns=['pos', 'gap_size']).set_index('pos')[
        'gap_size']
    all_gap_ref_positions_max.to_csv(
        os.path.join(output_dir, region_atlas.replace('dereplication.fa', 'reference_gap_size_per_position.csv')))
    reference_seq = list(SeqIO.parse(ref_file, 'fasta'))[0]
    ref_nuc = list(reference_seq.seq)
    seq_id = reference_seq.id
    seq_description = reference_seq.description
    ref_with_gaps = pd.Series(['N'] * (len(ref_nuc) + all_gap_ref_positions_max.sum()))
    correction_gap_pos = {}
    for gap_n, gap_pos in enumerate(all_gap_ref_positions_max.index):
        correction_gap_pos[gap_pos] = gap_pos + all_gap_ref_positions_max[:gap_n].sum().astype(int)
        ref_with_gaps[
        correction_gap_pos[gap_pos]:(correction_gap_pos[gap_pos] + all_gap_ref_positions_max[gap_pos])] = '-'
    ref_with_gaps[ref_with_gaps == 'N'] = ref_nuc
    all_aligned_records.append(
        SeqRecord.SeqRecord(Seq.Seq("".join(list(ref_with_gaps.values))), id=seq_id, description=seq_description))
    double_gap_reference_n = 0
    for variant_file in other_variants:
        aligned_seqs = list(SeqIO.parse(variant_file + '.needle', 'fasta'))
        currect_ref_record_seq = pd.Series(list(aligned_seqs[0].seq))
        currect_ref_record_gap = list(currect_ref_record_seq[currect_ref_record_seq == '-'].index.values)
        for gap_n, gap_pos in enumerate(currect_ref_record_gap):
            currect_ref_record_gap[gap_n] = gap_pos - gap_n  # gaps make the sequence longer which we correct
        how_many_times = pd.Series(currect_ref_record_gap).value_counts().astype(int)
        need_to_put_gap = []
        aligned_with_gaps = pd.Series(['N'] * len(ref_with_gaps))
        for gap_pos in correction_gap_pos.keys():
            if gap_pos not in currect_ref_record_gap:
                # need_to_put_gap.append(correction_gap_pos[gap_pos])
                how_map_gaps_in_pos = all_gap_ref_positions_max[gap_pos]
                corrected_pos = correction_gap_pos[gap_pos]
                aligned_with_gaps[corrected_pos:corrected_pos + how_map_gaps_in_pos] = '-'
            else:
                how_map_gaps_in_pos = all_gap_ref_positions_max[gap_pos] - how_many_times[gap_pos]
                assert how_map_gaps_in_pos >= 0
                if how_map_gaps_in_pos > 0:
                    corrected_pos = correction_gap_pos[gap_pos]
                    aligned_with_gaps[corrected_pos:corrected_pos + how_map_gaps_in_pos] = '-'
        aligned_nuc = list(aligned_seqs[1].seq)
        aligned_id = aligned_seqs[1].id
        aligned_description = aligned_seqs[1].description
        try:
            aligned_with_gaps[aligned_with_gaps == 'N'] = aligned_nuc
            all_aligned_records.append(
                SeqRecord.SeqRecord(Seq.Seq("".join(list(aligned_with_gaps.values))), id=aligned_id,
                                    description=aligned_description))
            if len(all_aligned_records)%5000==0:
                print('Wrote 5000 to %s'%output_region_alignment)
                with open(output_region_alignment, 'w') as f_h:
                    SeqIO.write(all_aligned_records, f_h, 'fasta')
        except ValueError:
            print(variant_file)
            raise ValueError
    with open(output_region_alignment, 'w') as f_h:
        SeqIO.write(all_aligned_records, f_h, 'fasta')
    print(os.path.basename(region_atlas), 'missing %s seqs' % double_gap_reference_n)

def run_needle_on_all_variants(output_dir,exe_path=''):#/oak/stanford/groups/pritch/users/daphna/tools/EMBOSS-6.6.0/emboss/'):
    for region_atlas in glob.glob(os.path.join(output_dir, '*.dereplication.fa')):
        tmp_name=os.path.basename(region_atlas).replace('.dereplication.fa','_tmp')
        tmp_dir=os.path.join(os.path.dirname(region_atlas),tmp_name)
        if not os.path.exists(tmp_dir):
            os.mkdir(os.path.join(os.path.dirname(region_atlas),tmp_name))
        run_needle_on_region_variants(region_atlas, output_dir, tmp_dir, exe_path=exe_path)

def msa_local():
    for region_atlas in glob.glob(os.path.join(output_dir, '*.dereplication.fa')):
        region = os.path.basename(region_atlas).replace('.dereplication.fa', '')
        msa_output = os.path.join(output_dir, region + '.clustalo.fa')
        command='clustalo --in %s --out %s'%(region_atlas,msa_output)
        os.system(command)

def get_atlas_needle_after_MSA(atlas_path,msa_path,output_dir):
    for atlas_region_f in glob.glob(os.path.join(atlas_path,'*atlas.fa')):
        atlas_ids = []
        if os.path.basename(atlas_region_f)=='all_region.atlas.fa':
            continue
        for record in SeqIO.parse(atlas_region_f,'fasta'):
            new_id = record.id
            new_id = new_id.split('_size_')[0]
            if '|' in new_id:
                new_id=new_id.split('|')[1]
            atlas_ids.append(new_id)
        region = os.path.basename(atlas_region_f).split('.atlas.')[0]
        found_records=[]
        for record in SeqIO.parse(os.path.join(msa_path,'all.MSA.%s.fa'%region),'fasta'):
            if record.id in atlas_ids:
                found_records.append(record)
        assert len(found_records)==len(atlas_ids)
        with open(os.path.join(output_dir,'%s.Needle_SA.fa'%region),'w') as f_h:
            SeqIO.write(found_records,f_h,'fasta')

def make_intermediante_vcf(output_dir):
    def iter(seq):#remove all "-" but if empty put "-"
        new_seq=seq.replace('-', '')
        if new_seq=='':
            return '-'
        return new_seq
    res = []
    for region_atlas in glob.glob(os.path.join(output_dir, '*.Needle_SA.fa')):
        region = os.path.basename(region_atlas).replace('.Needle_SA.fa', '')
        region_records=list(SeqIO.parse(region_atlas,'fasta'))
        reference=str(region_records[0].seq)
        alternative_records=region_records[1:]
        start_pos=0
        end_pos=0
        reference_non_gap_pos=0
        while end_pos<len(reference):
            cur_nuc=reference[end_pos]
            while cur_nuc=='-' and (end_pos<len(reference)-1):
                end_pos+=1
                cur_nuc = reference[end_pos]
            end_pos += 1
            ref_seq="".join(list(reference)[start_pos:end_pos]).upper()
            for record in alternative_records:
                alt_seq="".join(list(str(record.seq))[start_pos:end_pos])
                if alt_seq!=ref_seq:
                    if reference_non_gap_pos==len(reference.replace('-','')): #End gap
                        ref_seq = "".join(list(reference)[start_pos-1:end_pos]).upper()
                        alt_seq = "".join(list(str(record.seq))[start_pos-1:end_pos])
                        reference_non_gap_pos-=1
                    res.append(['s45',region,start_pos,end_pos,ref_seq,alt_seq,
                                reference_non_gap_pos,len(reference.replace('-',''))])
            start_pos=end_pos
            reference_non_gap_pos+=1
    res_df=pd.DataFrame(res,columns=['#CHROM','INFO','internal_start','internal_end','REF','ALT','POS','REF_LEN']).drop_duplicates()
    res_df=res_df[['#CHROM','POS','REF','ALT','INFO']]
    res_df['REF']=res_df['REF'].apply(lambda x: iter(x))
    res_df['ALT']=res_df['ALT'].apply(lambda x: iter(x))
    res_df=res_df.drop_duplicates()
    return res_df

def make_vcf(res_df,output_dir,throw_boundary=False,get_region=True):
    res_df_new=res_df.copy()
    helix_pos=pd.read_csv(helix_pos_table)
    es_pos = pd.read_csv(es_pos_table)
    region_45s_pos=pd.concat([helix_pos,es_pos]).drop_duplicates().set_index('Region')
    res_df_new['Region_POS']=res_df_new['POS'].values
    if get_region:
        res_df_new['POS']=res_df_new.apply(lambda x: x['POS']+region_45s_pos.loc[x['INFO'],'Start'],1)
    else:
        res_df_new['POS']=res_df_new['Region_POS']
    res=[]
    for pos_45s,pos_df in res_df_new.groupby('POS'):
        new_alt=",".join(list(set(pos_df['ALT'].values)))
        length = len(list(set(pos_df['ALT'].values)))
        res.append(['s45',pos_45s,pos_df['REF'].values[0],new_alt,pos_df['INFO'].values[0],length,pos_df['Region_POS'].values[0]])
    united_pos_df=pd.DataFrame(res,columns=['#CHROM','POS','REF','ALT','INFO','ALT_LEN','REGION_POS'])
    if throw_boundary:
        united_pos_df=united_pos_df.loc[~united_pos_df['POS'].isin(region_45s_pos['Start'])]
        united_pos_df = united_pos_df.loc[~united_pos_df['POS'].isin(region_45s_pos['End'])]
    united_pos_df.to_csv(os.path.join(output_dir,'SNV.atlas.Needle_SA.vcf'),index=False,sep='\t')
    return united_pos_df

def calc_idels(obj):
    indels=obj['ALT'].apply(
        lambda x: 1 if x == '-' else pd.Series(x.split(',')).apply(lambda x: len(x) > 1 if x != '-' else True).astype(
            int).sum()).sum()
    variants=obj['ALT_LEN'].sum()
    return indels/variants


if __name__=='__main__':
    region=os.path.join(output_tmp_dir,sys.argv[1]+'.dereplication.fa')
    region_tmp_dir = os.path.join(output_tmp_dir, sys.argv[1]+'_tmp')
    if not os.path.exists(region_tmp_dir):
        os.mkdir(region_tmp_dir)
    run_needle_on_region_variants(region,output_tmp_dir,region_tmp_dir,'EMBOSS-6.6.0/emboss/')
