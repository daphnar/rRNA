import pandas as pd
import pysam
import sys

lookup_path='all_es.atlas.position.variants.csv'
atlas_path = 'atlas_%s.hESC_mono_polysome.needle_gap10_extenstion_05.with_raw_count.csv'

def count_matches(f_in,f_out,binary=True):
    res = []
    if binary:
        mapiter = pysam.AlignmentFile(f_in, 'r')
    else:
        mapiter = pysam.AlignmentFile(f_in, 'rb')
    for read in mapiter.fetch():
        if read.reference_name is not None: #Removing unmapped
            if read.get_tag('NM') == 0:
                res.append(read.reference_name)
    pd.Series(res).value_counts().to_csv(f_out)

def get_nuc_freq(input_es_mapped_res,output_nuc_freq_mapped_res):
    mapped_input = pd.read_csv(input_es_mapped_res,index_col=0)
    mapped_input=mapped_input.rename(columns={'0':'count'})
    mapped_input['region']=mapped_input.index.map(lambda x: x.split(':')[1])
    found_variants = []
    translation = pd.read_csv(lookup_path).set_index('ID')
    for region,region_df in mapped_input.groupby('region'):
        region_ref=None
        all_region_positions=[]
        for id in region_df.index:
            list_of_positions = eval(translation.loc[id,'positions'])
            for idx,position in enumerate(list_of_positions):
                all_region_positions.append(position)
            if len(list_of_positions)==0:
                region_ref = id
        if region_ref is not None:
            for position in list_of_positions:
                found_variants.append([region, position, 0, region_df.loc[region_ref, 'count']])
        for id in region_df.index:
            list_of_positions = eval(translation.loc[id,'positions'])
            list_of_digits = eval(translation.loc[id,'variants'])
            for idx,position in enumerate(list_of_positions):
                found_variants.append([region,position,list_of_digits[idx],region_df.loc[id,'count']])
            #THIS COUNTS TWICE THE REFERENCE!
            for ref_not_count_positions in set(all_region_positions).difference(set(list_of_positions)):
                found_variants.append([region,ref_not_count_positions,0,region_df.loc[id,'count']])

    res_digits = pd.DataFrame(found_variants,columns=['region','position','digit','read_count'])
    res_digits_total=res_digits.groupby(['region','position','digit']).sum()
    res_digits_total=res_digits_total.reset_index()
    #res_digits_total['gene']=res_digits_total['region'].apply(lambda x: '18s' if x.startswith('SSU') else '28s')
    res_digits_total['gene'] = res_digits_total['region'].apply(lambda x: '18s' if x[-1] == 's' else '28s')

    res_nucs=[]
    for subunit,subunit_df in res_digits_total.groupby('gene'):
        translate_digit_to_nuc = pd.read_csv(atlas_path%subunit).set_index(['position', 'digit'])
        subunit_df = subunit_df.reset_index().set_index(['gene','position','digit'])
        for variant in subunit_df.index:
            position=variant[1]
            digit=variant[2]
            try:
                nuc=translate_digit_to_nuc.loc[(position,digit)]['variant']
                reference=translate_digit_to_nuc.loc[(position,digit)]['reference']
                if nuc==reference:
                    nuc_type='Reference'
                elif (len(nuc) == 1) and (nuc != '-'):
                    nuc_type = 'SNV'
                else:
                    nuc_type = 'indel'
                res_nucs.append([subunit, subunit_df.loc[variant, 'region'], position, nuc, digit,nuc_type,
                                 subunit_df.loc[variant, 'read_count']])
            except TypeError: #If a variant is not in the atlas but is found in the data
                pass

    res_nucs = pd.DataFrame(res_nucs, columns=['gene','region', 'position', 'nuc','digit','variant','read_count'])
    region_norm_factor=res_nucs.groupby(['region','position'])['read_count'].sum()
    norm_variant_freq=[]
    for idx in res_nucs.index:
        norm_variant_freq.append(res_nucs.loc[idx, 'read_count'] / region_norm_factor.loc[(res_nucs.loc[idx, 'region'],
                                                                                           res_nucs.loc[idx, 'position'])])

    res_nucs['region_freq']= norm_variant_freq
    res_nucs=res_nucs.drop(['region','digit'],1)
    res_nucs.to_csv(output_nuc_freq_mapped_res)

if __name__=='__main__':
    f_in = sys.argv[1]
    f_out='mapped_result_atlas_expand100_Helix.csv'
    binary=False
    count_matches(f_in, f_out, binary)
    f_nuc_freq = sys.argv[2]#'mapped_nucleotide_result_atlas.csv'
    get_nuc_freq(f_out, f_nuc_freq)

