#!/usr/bin/env python3
"""Import the squamata (non-snake) multi-genome SINEderella run into SINE-KB.

This script:
1. Loads the squamata sinekb_report.json (79 genomes, 9 subfamilies)
2. Imports per-genome copy counts using models.import_sinekb_report()
3. Extends the taxonomy tree with non-snake squamata species
4. Rebuilds the consensus bank

Usage:
    python import_squamata_run.py
"""

import json
import sys
from pathlib import Path

# Add sine-kb to path
sys.path.insert(0, str(Path(__file__).parent / "sine-kb"))
import models

# ═══════════════════════════════════════════════════════════════════════════════
# Taxonomy for 79 non-snake squamata species
# ═══════════════════════════════════════════════════════════════════════════════

SQUAMATA_TAXONOMY = {
    # Amphisbaenia
    "Rhineura_floridana": ("Rhineuridae", "Amphisbaenia"),
    # Anguimorpha
    "Dopasia_gracilis": ("Anguidae", "Anguimorpha"),
    "Elgaria_multicarinata": ("Anguidae", "Anguimorpha"),
    "Anniella_stebbinsi": ("Anniellidae", "Anguimorpha"),
    "Heloderma_charlesbogerti": ("Helodermatidae", "Anguimorpha"),
    "Lanthanotus_borneensis": ("Lanthanotidae", "Anguimorpha"),
    "Shinisaurus_crocodilurus": ("Shinisauridae", "Anguimorpha"),
    "Varanus_komodoensis": ("Varanidae", "Anguimorpha"),
    "Varanus_acanthurus": ("Varanidae", "Anguimorpha"),
    # Dibamia
    "Dibamus_smithi": ("Dibamidae", "Dibamia"),
    # Gekkota - Carphodactylidae
    "Nephrurus_levis": ("Carphodactylidae", "Gekkota"),
    # Gekkota - Diplodactylidae
    "Correlophus_ciliatus": ("Diplodactylidae", "Gekkota"),
    "Correlophus_sarasinorum": ("Diplodactylidae", "Gekkota"),
    # Gekkota - Eublepharidae
    "Coleonyx_brevis": ("Eublepharidae", "Gekkota"),
    "Coleonyx_elegans": ("Eublepharidae", "Gekkota"),
    "Eublepharis_macularius": ("Eublepharidae", "Gekkota"),
    # Gekkota - Gekkonidae
    "Christinus_marmoratus": ("Gekkonidae", "Gekkota"),
    "Gekko_gecko": ("Gekkonidae", "Gekkota"),
    "Gekko_japonicus": ("Gekkonidae", "Gekkota"),
    "Gehyra_mutilata": ("Gekkonidae", "Gekkota"),
    "Heteronotia_binoei": ("Gekkonidae", "Gekkota"),
    "Hemidactylus_brookii": ("Gekkonidae", "Gekkota"),
    "Hemidactylus_flaviviridis": ("Gekkonidae", "Gekkota"),
    "Hemidactylus_frenatus": ("Gekkonidae", "Gekkota"),
    "Hemidactylus_mabouia": ("Gekkonidae", "Gekkota"),
    "Hemidactylus_macropholis": ("Gekkonidae", "Gekkota"),
    "Lepidodactylus_listeri": ("Gekkonidae", "Gekkota"),
    "Lepidodactylus_lugubris": ("Gekkonidae", "Gekkota"),
    "Mediodactylus_kotschyi": ("Gekkonidae", "Gekkota"),
    "Phelsuma_laticauda": ("Gekkonidae", "Gekkota"),
    "Paroedura_picta": ("Gekkonidae", "Gekkota"),
    "Paroedura_stumpffi": ("Gekkonidae", "Gekkota"),
    # Gekkota - Phyllodactylidae
    "Asaccus_caudivolvulus": ("Phyllodactylidae", "Gekkota"),
    "Phyllodactylus_wirshingi": ("Phyllodactylidae", "Gekkota"),
    "Tarentola_annularis": ("Phyllodactylidae", "Gekkota"),
    "Tarentola_mauritanica": ("Phyllodactylidae", "Gekkota"),
    "Thecadactylus_rapicauda": ("Phyllodactylidae", "Gekkota"),
    # Gekkota - Pygopodidae
    "Lialis_burtonis": ("Pygopodidae", "Gekkota"),
    "Pygopus_nigriceps": ("Pygopodidae", "Gekkota"),
    # Gekkota - Sphaerodactylidae
    "Aristelliger_praesignis": ("Sphaerodactylidae", "Gekkota"),
    "Euleptes_europaea": ("Sphaerodactylidae", "Gekkota"),
    "Gonatodes_alexandermendesi": ("Sphaerodactylidae", "Gekkota"),
    "Sphaerodactylus_glaucus": ("Sphaerodactylidae", "Gekkota"),
    "Sphaerodactylus_inigoi": ("Sphaerodactylidae", "Gekkota"),
    "Sphaerodactylus_townsendi": ("Sphaerodactylidae", "Gekkota"),
    "Teratoscincus_roborowskii": ("Sphaerodactylidae", "Gekkota"),
    # Iguania
    "Laudakia_wui": ("Agamidae", "Iguania"),
    "Pogona_vitticeps": ("Agamidae", "Iguania"),
    "Chamaeleo_calyptratus": ("Chamaeleonidae", "Iguania"),
    "Furcifer_pardalis": ("Chamaeleonidae", "Iguania"),
    "Basiliscus_vittatus": ("Corytophanidae", "Iguania"),
    "Gambelia_wislizenii": ("Crotaphytidae", "Iguania"),
    "Anolis_carolinensis": ("Dactyloidae", "Iguania"),
    "Cyclura_pinguis": ("Iguanidae", "Iguania"),
    "Iguana_iguana": ("Iguanidae", "Iguania"),
    "Brachylophus_fasciatus": ("Iguanidae", "Iguania"),
    # Laterata
    "Gallotia_galloti": ("Lacertidae", "Laterata"),
    "Latastia_doriai": ("Lacertidae", "Laterata"),
    "Podarcis_lilfordi": ("Lacertidae", "Laterata"),
    "Podarcis_cretensis": ("Lacertidae", "Laterata"),
    "Podarcis_muralis": ("Lacertidae", "Laterata"),
    "Lacerta_agilis": ("Lacertidae", "Laterata"),
    "Zootoca_vivipara": ("Lacertidae", "Laterata"),
    "Tretioscincus_oriximinensis": ("Gymnophthalmidae", "Laterata"),
    "Aspidoscelis_tigris": ("Teiidae", "Laterata"),
    "Salvator_merianae": ("Teiidae", "Laterata"),
    # Scinciformata
    "Cryptoblepharus_egeriae": ("Scincidae", "Scinciformata"),
    "Hemicordylus_capensis": ("Cordylidae", "Scinciformata"),
    "Lerista_edwardsae": ("Scincidae", "Scinciformata"),
    "Marisora_unimarginata": ("Scincidae", "Scinciformata"),
    "Spondylurus_nitidus": ("Scincidae", "Scinciformata"),
    "Tiliqua_scincoides": ("Scincidae", "Scinciformata"),
    "Xantusia_vigilis": ("Xantusiidae", "Scinciformata"),
}

# Duplicate assemblies in the multi-run (map variant -> canonical)
DUPLICATES = {
    "Correlophus_ciliatus_local": "Correlophus_ciliatus",
    "Eublepharis_macularius_local": "Eublepharis_macularius",
    "Gekko_japonicus_local": "Gekko_japonicus",
    "Gehyra_mutilata_hifi": "Gehyra_mutilata",
    "Gehyra_mutilata_local": "Gehyra_mutilata",
    "Heteronotia_binoei_local": "Heteronotia_binoei",
    "Thecadactylus_rapicauda_local": "Thecadactylus_rapicauda",
}

# ═══════════════════════════════════════════════════════════════════════════════
# Build the full taxonomy tree (Serpentes + non-snake Squamata)
# Kept from import_snake_run.py and extended with outgroups
# ═══════════════════════════════════════════════════════════════════════════════

def build_full_tree():
    """Build the complete Squamata taxonomy tree."""
    # ── Snake portion (from import_snake_run.py) ──
    serpentes = {
        "name": "Serpentes",
        "children": [
            {"name": "Typhlopidae", "children": [
                {"taxon": "Anilios_bituberculatus", "common_name": "Southern blind snake"},
                {"taxon": "Indotyphlops_braminus", "common_name": "Brahminy blind snake"},
            ]},
            {"name": "Aniliidae", "children": [
                {"taxon": "Anilius_scytale", "common_name": "Red pipe snake"},
            ]},
            {"name": "Uropeltidae", "children": [
                {"taxon": "Rhinophis_saffragamus", "common_name": "Sri Lanka earth snake"},
            ]},
            {"name": "Cylindrophiidae", "children": [
                {"taxon": "Cylindrophis_ruffus", "common_name": "Red-tailed pipe snake"},
            ]},
            {"name": "Xenopeltidae", "children": [
                {"taxon": "Xenopeltis_unicolor", "common_name": "Sunbeam snake"},
            ]},
            {"name": "Loxocemidae", "children": [
                {"taxon": "Loxocemus_bicolor", "common_name": "Mexican burrowing python"},
            ]},
            {"name": "Boidae", "children": [
                {"name": "Erycinae", "children": [
                    {"taxon": "Eryx_jayakari", "common_name": "Arabian sand boa"},
                    {"taxon": "Eryx_tataricus", "common_name": "Tatar sand boa"},
                ]},
                {"name": "Boinae", "children": [
                    {"taxon": "Boa_constrictor", "common_name": "Boa constrictor"},
                    {"taxon": "Corallus_caninus", "common_name": "Emerald tree boa"},
                ]},
                {"name": "Candoiinae", "children": [
                    {"taxon": "Candoia_aspera", "common_name": "Viper boa"},
                ]},
                {"name": "Charininae", "children": [
                    {"taxon": "Charina_bottae", "common_name": "Rubber boa"},
                    {"taxon": "Lichanura_trivirgata", "common_name": "Rosy boa"},
                ]},
            ]},
            {"name": "Pythonidae", "children": [
                {"taxon": "Python_bivittatus", "common_name": "Burmese python"},
                {"taxon": "Malayopython_reticulatus", "common_name": "Reticulated python"},
                {"taxon": "Liasis_olivaceus", "common_name": "Olive python"},
                {"taxon": "Morelia_carinata", "common_name": "Rough-scaled python"},
                {"taxon": "Morelia_viridis", "common_name": "Green tree python"},
                {"taxon": "Simalia_boeleni", "common_name": "Boelen's python"},
                {"taxon": "Simalia_tracyae", "common_name": "Tracy's scrub python"},
            ]},
            {"name": "Acrochordidae", "children": [
                {"taxon": "Acrochordus_granulatus", "common_name": "Little file snake"},
            ]},
            {"name": "Caenophidia", "children": [
                {"name": "Colubridae", "children": [
                    {"taxon": "Arizona_elegans", "common_name": "Glossy snake"},
                    {"taxon": "Elaphe_schrenckii", "common_name": "Amur rat snake"},
                ]},
                {"name": "Elapidae", "children": [
                    {"taxon": "Naja_naja", "common_name": "Indian cobra"},
                    {"taxon": "Hydrophis_major", "common_name": "Olive-headed sea snake"},
                ]},
                {"name": "Homalopsidae", "children": [
                    {"taxon": "Hypsiscopus_plumbea", "common_name": "Plumbeous water snake"},
                    {"taxon": "Myanophis_thanlyinensis", "common_name": "Thanhlyin water snake"},
                ]},
                {"name": "Viperidae", "children": [
                    {"taxon": "Azemiops_feae", "common_name": "Fea's viper"},
                    {"taxon": "Crotalus_horridus", "common_name": "Timber rattlesnake"},
                    {"taxon": "Vipera_latastei", "common_name": "Lataste's viper"},
                ]},
                {"name": "Xenodermidae", "children": [
                    {"taxon": "Achalinus_rufescens", "common_name": "Burrowing odd-scaled snake"},
                ]},
            ]},
        ],
    }

    # ── Non-snake Squamata (outgroups) ──
    # Group by clade → family → species
    from collections import OrderedDict
    clades = OrderedDict()
    for species, (fam, clade) in sorted(SQUAMATA_TAXONOMY.items(), key=lambda x: (x[1][1], x[1][0], x[0])):
        if clade not in clades:
            clades[clade] = OrderedDict()
        if fam not in clades[clade]:
            clades[clade][fam] = []
        clades[clade][fam].append(species)

    outgroup_children = []
    for clade, families in clades.items():
        clade_children = []
        for fam, spp in families.items():
            fam_children = [{"taxon": sp} for sp in spp]
            clade_children.append({"name": fam, "children": fam_children})
        outgroup_children.append({"name": clade, "children": clade_children})

    tree = [
        {
            "name": "Reptilia",
            "children": [
                {
                    "name": "Squamata",
                    "children": [serpentes] + outgroup_children,
                }
            ],
        }
    ]
    return tree


def main():
    # Initialize sine-kb data layer
    data_dir = str(Path(__file__).parent / "sine-kb" / "data")
    models.init(data_dir)

    report_path = Path(__file__).parent / "server_import" / "squamata_sinekb_report.json"
    report = json.load(open(str(report_path)))

    print("=== Importing squamata (non-snake) multi-genome run ===")
    print("Report: {} genomes, {} subfamilies".format(len(report['genomes']), len(report['subfamilies'])))
    print()

    # ── Pre-import stats ─────────────────────────────────────────────────
    st = models.stats()
    print("Before import:")
    print("  Families:    {}".format(st['families']))
    print("  Subfamilies: {}".format(st['subfamilies']))
    print("  Instances:   {}".format(st['instances']))
    print("  Taxa:        {}".format(st['taxa']))
    print()

    # ── Import using built-in import_sinekb_report ───────────────────────
    summary = models.import_sinekb_report(report)

    print("Import results:")
    print("  Subfamilies created: {}".format(summary['subfamilies_created']))
    print("  Subfamilies updated: {}".format(summary['subfamilies_updated']))
    print("  Instances created:   {}".format(summary['instances_created']))
    print("  Instances updated:   {}".format(summary['instances_updated']))
    if summary.get('skipped'):
        print("  Skipped: {}".format(summary['skipped']))
    print()

    # ── Mark duplicate assemblies in instance notes ──────────────────────
    print("Marking duplicate assemblies...")
    dup_count = 0
    for variant, canonical in DUPLICATES.items():
        tslug_var = models.slug(variant)
        tslug_can = models.slug(canonical)
        # Find all instances for the variant
        fdir = Path(data_dir) / "families"
        for fam in sorted(fdir.iterdir()):
            if not fam.is_dir():
                continue
            sbase = fam / "subfamilies"
            if not sbase.exists():
                continue
            for sf in sorted(sbase.iterdir()):
                ip = sf / "taxa" / tslug_var / "instance.yaml"
                if ip.exists():
                    inst = models.get_instance(fam.name, sf.name, tslug_var)
                    if inst and not inst.get("notes"):
                        inst["notes"] = "Alternate assembly for {}".format(canonical)
                        models.save_instance(fam.name, sf.name, tslug_var, inst)
                        dup_count += 1
    print("  {} instances annotated as duplicates".format(dup_count))
    print()

    # ── Update taxonomy tree ─────────────────────────────────────────────
    full_tree = build_full_tree()
    models.save_tree(full_tree)
    leaves = models.get_leaf_taxa(full_tree)
    print("Taxonomy tree updated: {} species".format(len(leaves)))
    print()

    # ── Rebuild consensus bank ───────────────────────────────────────────
    n = models.build_bank()
    print("Consensus bank rebuilt: {} sequences".format(n))
    print()

    # ── Post-import stats ────────────────────────────────────────────────
    st = models.stats()
    print("=== SINE-KB Stats (after import) ===")
    print("  Families:    {}".format(st['families']))
    print("  Subfamilies: {}".format(st['subfamilies']))
    print("  Instances:   {}".format(st['instances']))
    print("  Taxa:        {}".format(st['taxa']))
    print()
    print("Done!")


if __name__ == "__main__":
    main()
