import hail as hl
import argparse

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
                          .when(((mt.alleles[0] == ss_loc.ref) &
                                 (mt.alleles[1] == ss_loc.alt)) |
                                ((flip_text(mt.alleles[0]) == ss_loc.ref) &
                                 (flip_text(mt.alleles[1]) == ss_loc.alt)),
                                (-1 * ss_loc.beta))
                          .when(((mt.alleles[0] == ss_loc.alt) &
                                 (mt.alleles[1] == ss_loc.ref)) |
                                ((flip_text(mt.alleles[0]) == ss_loc.alt) &
                                 (flip_text(mt.alleles[1]) == ss_loc.ref)),
                                ss_loc.beta)
                          .or_missing()}
                          )
    return(mt)

def main(args):
    betas = ['beta_01', 'beta_1', 'beta_10', 'beta_100']
    spike_slab = hl.import_table('gs://armartin/mama/spike_slab/BBJ_UKB_hm3.chr22.cm.beta.txt', impute=True)
    spike_slab = spike_slab.key_by(**hl.parse_variant(spike_slab.v))
    if args.compute_true_phenotypes:
        # get the white british subset
        eur = hl.import_table('gs://phenotype_31063/ukb31063.gwas_covariates.both_sexes.tsv').key_by('s')

        # read in imputed data, subset to chr22
        mt = hl.read_matrix_table('gs://phenotype_31063/hail/imputed/ukb31063.dosage.autosomes.mt')
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval('22')])

        # annotate and filter imputed data to all sites with causal effects
        mt = mt.annotate_rows(ss=spike_slab[mt.row_key])
        mt = mt.filter_rows(hl.is_defined(mt.ss))

        # compute true PRS (i.e. phenotypes)
        annot_expr = {
            i: hl.agg.sum(mt.ss[i] * mt.dosage) for i in betas
        }

        # write out phenos for white British unrelated subset
        mt = mt.annotate_cols(**annot_expr)
        mt = mt.filter_cols(hl.is_defined(eur[mt.s]))
        mt.cols().write('gs://armartin/mama/spike_slab/BBJ_UKB_hm3.chr22.cm.beta.true_PRS.ht', stage_locally=True, overwrite=True)

    if args.run_gwas:
        # read back in PRS (now true phenotypes)
        phenos = hl.read_table('gs://armartin/mama/spike_slab/BBJ_UKB_hm3.chr22.cm.beta.true_PRS.ht').key_by('s')
        phenos.show()
        covariates = hl.import_table('gs://phenotype_31063/ukb31063.gwas_covariates.both_sexes.tsv', impute=True,
                                     types={'s': hl.tstr}).key_by('s')
        full_mt = hl.read_matrix_table('gs://phenotype_31063/hail/imputed/ukb31063.dosage.autosomes.mt')
        full_mt = full_mt.annotate_cols(**covariates[full_mt.s])
        full_mt = hl.filter_intervals(full_mt, [hl.parse_locus_interval('22')])

        # annotate and filter imputed data to all sites with causal effects
        full_mt = full_mt.annotate_rows(ss=spike_slab[full_mt.row_key])
        full_mt = full_mt.filter_rows(hl.is_defined(full_mt.ss))

        # subset to white British subset, get 10 sets of 10k and run a gwas for each of these w/ PCs as covs
        for i in range(10):
            subset_pheno = phenos.annotate(r=hl.rand_unif(0, 1))
            subset_pheno = subset_pheno.order_by(subset_pheno.r).add_index('global_idx').key_by('s')
            subset_pheno = subset_pheno.filter(subset_pheno.global_idx < 10000)
            mt = full_mt.annotate_cols(**subset_pheno[full_mt.s])
            mt = mt.annotate_rows(maf=hl.agg.mean(mt.dosage)/2)
            result_ht = hl.linear_regression_rows(
                y=[mt[i] for i in betas],
                x=mt.dosage,
                covariates=[1] + [mt['PC' + str(i)] for i in range(1, 21)],
                pass_through=['rsid', 'maf'])

            subset_pheno.export('gs://armartin/mama/spike_slab/UKB_hm3.chr22.cm.beta.true_PRS.gwas_inds_' + str(i) + '.tsv.gz')
            result_ht.write('gs://armartin/mama/spike_slab/UKB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + str(i) + '.ht', overwrite=True)

    if args.write_gwas:
        for i in range(10):
            result_ht = hl.read_table('gs://armartin/mama/spike_slab/UKB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + str(i) + '.ht')
            result_ht = result_ht.key_by()
            get_expr = {
                field + '_' + x: result_ht[field][i] for i, x in enumerate(betas) for field in ['beta', 'standard_error', 'p_value']
            }
            result_ht.select(chr=result_ht.locus.contig, pos=result_ht.locus.position, rsid=result_ht.rsid, ref=result_ht.alleles[0],
                             alt=result_ht.alleles[1], maf=result_ht.maf, n=result_ht.n, **get_expr)\
                .export('gs://armartin/mama/spike_slab/UKB_hm3.chr22.cm.beta.true_PRS.gwas_sumstat_' + str(i) + '.tsv.gz')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compute_true_phenotypes', action='store_true')
    parser.add_argument('--run_gwas', action='store_true')
    parser.add_argument('--write_gwas', action='store_true')
    args = parser.parse_args()
    main(args)

