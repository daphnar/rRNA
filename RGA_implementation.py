def run_needle_on_region_variants(region_atlas,output_dir,tmp_dir,exe_path=''):
    executable = exe_path + 'needle'
    output_region_alignment = region_atlas.replace('dereplication', 'Needle_SA')
    all_aligned_records = []
    seq_files = split_allseqs_to_single_files(region_atlas, tmp_dir)
    ref_file = seq_files[0]
    other_variants = seq_files[1:]
    all_gap_ref_positions = []
    for idx,variant_file in enumerate(other_variants):
        ##################OLD params:
        # command = '%s -asequence %s -bsequence %s -datafile EDNAFULL -gapopen 2' \
        #           ' -gapextend 10 -snucleotide1' \
        #           ' -snucleotide2 -outfile %s -aformat fasta' % (
        #           executable, ref_file, variant_file, variant_file + '.needle')
        ##################New params:
        command = '%s -asequence %s -bsequence %s -datafile EDNAFULL -gapopen 10' \
                  ' -gapextend 0.5 -snucleotide1' \
                  ' -snucleotide2 -outfile %s -aformat fasta' % (
                      executable, ref_file, variant_file, variant_file + '.needle')
        if not os.path.exists(variant_file + '.needle'):
            run_shell_command(command, False)
        if os.stat(variant_file + '.needle').st_size==0:
            run_shell_command(command, False)
        ref_record_s = pd.Series(list(list(SeqIO.parse(variant_file + '.needle', 'fasta'))[0].seq))
        ref_record_gap = list(ref_record_s[ref_record_s == '-'].index.values)
        for gap_n, gap_pos in enumerate(ref_record_gap):
            ref_record_gap[gap_n] = gap_pos - gap_n  # gaps make the sequence longer which we correct
        # all_gap_ref_positions+=ref_record_gap
        how_many_times = pd.Series(ref_record_gap).value_counts().astype(int)
        all_gap_ref_positions.append(how_many_times)
        if idx%5000==0:
            print('5000 done')
    # all_gap_ref_positions=list(set(all_gap_ref_positions))
    # all_gap_ref_positions.sort()
    # print('WTF %s'%(len(all_gap_ref_positions)))
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
    # #Todo: uncomment above and remove line below
    # all_gap_ref_positions_max = pd.read_csv(os.path.join(output_dir, region_atlas.replace('dereplication.fa', 'reference_gap_size_per_position.csv')),index_col=0)['gap_size']
    reference_seq = list(SeqIO.parse(ref_file, 'fasta'))[0]
    ref_nuc = list(reference_seq.seq)
    seq_id = reference_seq.id
    seq_description = reference_seq.description
    # ref_with_gaps=pd.Series(['N']*(len(ref_nuc)+len(all_gap_ref_positions)))
    ref_with_gaps = pd.Series(['N'] * (len(ref_nuc) + all_gap_ref_positions_max.sum()))
    correction_gap_pos = {}
    # for gap_n,gap_pos in enumerate(all_gap_ref_positions):
    for gap_n, gap_pos in enumerate(all_gap_ref_positions_max.index):
        # correction_gap_pos[gap_pos]=gap_pos+gap_n #gaps make the sequence longer which we correct
        correction_gap_pos[gap_pos] = gap_pos + all_gap_ref_positions_max[:gap_n].sum().astype(int)
        ref_with_gaps[
        correction_gap_pos[gap_pos]:(correction_gap_pos[gap_pos] + all_gap_ref_positions_max[gap_pos])] = '-'
    # ref_with_gaps[correction_gap_pos.values()] = '-'
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
        #good=True#REMOVE
        for gap_pos in correction_gap_pos.keys():
            if gap_pos not in currect_ref_record_gap:
                # need_to_put_gap.append(correction_gap_pos[gap_pos])
                how_map_gaps_in_pos = all_gap_ref_positions_max[gap_pos]
                corrected_pos = correction_gap_pos[gap_pos]
                aligned_with_gaps[corrected_pos:corrected_pos + how_map_gaps_in_pos] = '-'
            else:
                how_map_gaps_in_pos = all_gap_ref_positions_max[gap_pos] - how_many_times[gap_pos]
                assert how_map_gaps_in_pos >= 0
                # if how_map_gaps_in_pos < 0:
                #     good=False
                #     break
                if how_map_gaps_in_pos > 0:
                    corrected_pos = correction_gap_pos[gap_pos]
                    aligned_with_gaps[corrected_pos:corrected_pos + how_map_gaps_in_pos] = '-'
        # if not good:
        #     continue
        # aligned_with_gaps[need_to_put_gap]='-'
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
