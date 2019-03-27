import hail as hl
import pickle
import time
import argparse
from pprint import pprint

def flip_text(base):
    """
    :param StringExpression base: Expression of a single base
    :return: StringExpression of flipped base
    :rtype: StringExpression
    """
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('C', 'G')
            .when('G', 'C')
            .default(base))


def annotate_beta(mt, ss_loc):
    mt = mt.annotate_rows(**{
        'beta': hl.case()
                          .when(((mt.alleles[0] == ss_loc.A1) &
                                 (mt.alleles[1] == ss_loc.A2)) |
                                ((flip_text(mt.alleles[0]) == ss_loc.A1) &
                                 (flip_text(mt.alleles[1]) == ss_loc.A2)),
                                (-1 * ss_loc.MAMA_BETA))
                          .when(((mt.alleles[0] == ss_loc.A2) &
                                 (mt.alleles[1] == ss_loc.A1)) |
                                ((flip_text(mt.alleles[0]) == ss_loc.A2) &
                                 (flip_text(mt.alleles[1]) == ss_loc.A1)),
                                ss_loc.MAMA_BETA)
                          .or_missing()}
                          )
    return(mt)


def specific_clumps(filename):
    clump = hl.import_table(filename, delimiter='\s+', min_partitions=10, types={'P': hl.tfloat})
    clump_dict = clump.aggregate(hl.dict(hl.agg.collect(
        (hl.locus(hl.str(clump.CHR), hl.int(clump.BP)),
        True)
    )), _localize=False)
    return clump_dict


def main(args):
    ########################################################################
    ### initialize
    which_beta = 'beta' + args.which_beta
    if args.method == 'metal':
        end_dir = 'metal/'
        clumps = args.dirname + end_dir + 'BBJ_UKBB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_beta_' + args.which_beta + '_' + args.iter + '.metal.clumped'
        ss_filename = args.dirname + end_dir + 'BBJ_UKBB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_beta_' + args.which_beta + '_' + args.iter + '.tsv'
        out_base = args.dirname + end_dir + which_beta + '_draw_' + args.iter + '_spike_' + args.which_beta + '_metal_PRS'
    elif args.method == 'mama':
        ld = args.ld + '/'
        analysis = args.analysis + '/'
        end_dir = 'mama/ld_true/'
        clumps = args.dirname + end_dir + ld + analysis + 'draw_' + args.iter + '_spike_' + args.which_beta + '_mama_2.clumped'
        ss_filename = args.dirname + end_dir + ld + analysis + 'draw_' + args.iter + '_spike_' + args.which_beta + '_mama_2.txt'
        # this ss_filename has different headers

        out_base = args.dirname + end_dir + ld + analysis + 'draw_' + args.iter + '_spike_' + args.which_beta + '_mama_2_PRS'
    else:
        end_dir = 'ukbb_only/'
        clumps = args.dirname + end_dir + 'UKB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + args.iter + '_beta' + args.which_beta + '.clumped'
        ss_filename = args.dirname + end_dir + 'UKB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + args.iter + '.tsv.gz'
        out_base = args.dirname + end_dir + 'UKB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + args.iter + '_beta' + args.which_beta + '_gwas_PRS'

    clump_table_location = args.dirname + 'keytables/ukb-' + args.basename + '-pt-sumstats-locus-allele-keyed.kt'

    contigs = {'0{}'.format(x):str(x) for x in range(1, 10)}

    bgen_files = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr22_v3.bgen'

    start = time.time()
    # large block size because we read very little data (due to filtering & ignoring genotypes)
    hl.init(branching_factor=10, min_block_size=2000)
    # set min_block_size only in import_bgen


    ################################################################################
    ### set up the sumstats table (chr, bp for union SNPs)
    if args.read_clumps:
        clump_file = hl.import_table(clumps,
                                     delimiter='\s+',
                                     impute=True)
        clump_file = clump_file.select(locus=hl.locus(hl.str(clump_file.CHR), clump_file.BP))
        clump_file = clump_file.key_by('locus')
        clump_file.write(clump_table_location, overwrite=True)

    clump_file = hl.read_table(clump_table_location)

    ################################################################################
    ### Write ss info, process so sumstats are uniform across MAMA, METAL, and gwas
    if args.ss_tables:
        #ss = hl.import_table(args.dirname + args.basename + '.tsv.gz',
        ss = hl.import_table(ss_filename,
                             #'BBJ_UKBB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + args.which_beta + 'beta_01_9.tsv'
                                                      #  # for mama case
                             delimiter='\s+',
                             impute=True,
                             types={'BP': hl.tint})
        if args.method != 'mama' and args.method != 'metal':
            ss = ss.rename({'chr': 'CHR', 'pos': 'BP', 'rsid': 'SNP', 'ref': 'A1', 'alt': 'A2',
                            'maf': 'FRQ', 'p_value_beta_' + args.which_beta: 'MAMA_PVAL',
                            'standard_error_beta_' + args.which_beta: 'MAMA_SE',
                            'beta_beta_' + args.which_beta: 'MAMA_BETA'})
        ss = ss.key_by(locus=hl.locus(hl.str(ss.CHR), hl.int(ss.BP))).repartition(200)
        ss = ss.annotate(A1=ss.A1.upper(), A2=ss.A2.upper())

        ss.write(args.dirname + args.basename + '_sep.ht', True)

    ss = hl.read_table(args.dirname + args.basename + '_sep.ht')

    ################################################################################
    ### Run the PRS using phenotype-specific clump variants
    if args.write_bgen:
        mt_all = hl.import_bgen(bgen_files,
                                ['dosage'],
                                sample_file='gs://phenotype_31063/ukb31063.autosomes.sample',
                                variants=clump_file.locus)

        samples = hl.import_table(args.dirname + 'ukb_not_in_simulation_rand5000.inds', types={'s': hl.tstr}).key_by('s')
        mt = mt_all.filter_cols(hl.is_defined(samples[mt_all.s]))

        mt.repartition(5000, shuffle=False).write(args.dirname + args.basename + '.mt', True)


    mt = hl.read_matrix_table(args.dirname + args.basename + '.mt')
    true_ss = hl.read_table(args.dirname + 'BBJ_UKB_hm3.chr22.cm.beta.true_PRS.ht')

    """
    To add:
    - Also fix nt1/nt2 to A1 and A2 (check) from sumstats.
    """
    # filter to only samples held out from GWAS
    mt = mt.annotate_rows(ss=ss[mt.locus])
    mt = annotate_beta(mt, mt.ss)

    p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5, 's10': 1}

    pheno_clump = specific_clumps(clumps)

    mt = mt.filter_rows(pheno_clump.get(mt.locus, False))
    print(mt.count())

    # divide by sd's of frequencies to get standardized betas back to allelic scale for MAMA betas (only, not METAL)
    # sqrt(2pq)
    if args.betas_are_standardized:
        annot_expr = {
            k: hl.agg.sum(mt.beta / hl.sqrt(2 * hl.float(mt.ss.FRQ) * (1-hl.float(mt.ss.FRQ))) * mt.dosage * hl.int(mt.ss.MAMA_PVAL < v))
            for k, v in p_max.items()}
    else:
        annot_expr = {
            k: hl.agg.sum(mt.beta * mt.dosage * hl.int(mt.ss.MAMA_PVAL < v))
            for k, v in p_max.items()}

    mt = mt.annotate_cols(**annot_expr, **true_ss[mt.s])

    mt.key_cols_by().cols().write(out_base + '.ht', stage_locally=True, overwrite=True)
    ht = hl.read_table(out_base + '.ht')

    output_location = out_base + '.txt.bgz'
    ht.export(output_location)
    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--read_clumps', action='store_true')
    parser.add_argument('--ss_tables', action='store_true')
    parser.add_argument('--write_bgen', action='store_true')
    parser.add_argument('--which_beta', default='01')
    parser.add_argument('--iter')
    parser.add_argument('--method')
    parser.add_argument('--ld')
    parser.add_argument('--analysis')
    parser.add_argument('--dirname', default='gs://armartin/disparities/bbj/') # gs://armartin/disparities/ukbb/
    parser.add_argument('--basename', default='BBJ_holdout_gwas_')  # pheno_31063_holdout_gwas_
    parser.add_argument('--betas_are_standardized', action='store_true')

    args = parser.parse_args()
    main(args)
    #try_slack('@armartin', main, args)
