"""
Certain features (such as lookups) require large files that are not included with this package.

This module governs the building and downloading of such assets.
"""
from filefetcher import AssetCLI, AssetManager


# For routine use, instantiate a manager. It will locate cached copies of asset files.
# FIXME: Temporary download URL; will eventually move this to an S3 bucket
manager = AssetManager('zorp', 'https://csg-assets.occsci.com/zorp/manifest.json')  # site hosts a manifest.json file


def set_recipes():
    """
    The Zorp use case assumes that all assets are built in advance, so we will only define recipes
        when explicitly asked to do so (such as in CLI mode)
    """

    from .loaders.make_rsid_lookup import MakeSnpToRsid

    B37_SAMPLE_GENES = (
        ('2', 21224301, 21266945),  # APOB
        ('19', 45409039, 45412650),  # APOE
        ('17', 41196312, 41277500),  # BRCA1
        ('13', 32889617, 32973809),  # BRCA2
        ('10', 135340300, 135352627),  # CYP2E1
        ('16', 53737875, 54148379),  # FTO
        ('15', 28356183, 28567313),  # HERC2
        ('1', 55505149, 55530526),  # PCSK9
        ('10', 114709978, 114927437),  # TCF7L2
    )

    B38_SAMPLE_GENES = (
        ('2', 21001429, 21044073),  # APOB
        ('19', 44905796, 44909395),  # APOE
        ('17', 43044295,  43125364),  # BRCA1
        ('13', 32315508, 32400268),  # BRCA2
        ('10', 133527363, 133539123),  # CYP2E1
        ('16', 53703963, 54121941),  # FTO
        ('15', 28111040, 28322179),  # HERC2
        ('1', 55039548, 55064853),  # PCSK9
        ('10', 112950247, 113167678),  # TCF7L2
    )

    manager.add_recipe(
        'snp_to_rsid',
        MakeSnpToRsid('GRCh37'),
        label='Find rsID information given chrom/pos/ref/alt',
        genome_build='GRCh37'
    )

    manager.add_recipe(
        'snp_to_rsid_test',
        MakeSnpToRsid('GRCh37', sample_regions=B37_SAMPLE_GENES),
        label='Find rsID information given chrom/pos/ref/alt',
        genome_build='GRCh37'
    )

    manager.add_recipe(
        'snp_to_rsid',
        MakeSnpToRsid('GRCh38'),
        label='Find rsID information given chrom/pos/ref/alt',
        genome_build='GRCh38'
    )

    manager.add_recipe(
        'snp_to_rsid_test',
        MakeSnpToRsid('GRCh38', sample_regions=B38_SAMPLE_GENES),
        label='Find rsID information given chrom/pos/ref/alt',
        genome_build='GRCh38'
    )


def main():
    """
    Main method. Register recipes and run the CLI.
    """
    set_recipes()

    cli = AssetCLI(manager)
    cli.run()


if __name__ == '__main__':
    main()
