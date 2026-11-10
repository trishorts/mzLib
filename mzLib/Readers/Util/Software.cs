namespace Readers
{
    public enum Software
    {
        Unspecified,
        MassSpecFile,

        // Deconvolution
        FLASHDeconv, // files tested were outputted from OpenMs3.0.0
        TopFD,       // files tested were outputted from v1.6.2
        IsoDec,

        // Search
        MetaMorpheus,
        MaxQuant,
        Toppic,
        MsFragger, // files tested were from fragpipe v21.1
        MsPathFinderT,
        Crux,
        Dinosaur,
        ProteomeDiscoverer,
        ProsightPD,
        Chimerys,
        DiaNn,     // files tested were from DIA-NN 1.8/1.9 long-format report.tsv
        Pytheas,
        Spectronaut, // files tested were from Spectronaut 10 (csv) and 19.2, and from MassIVE deposits of 2022-2025 (tsv, tab-text xls)

        // Quantification
        FlashLFQ,

        // Other
        Casanovo,
        MzIdentML, // a PSI standard written by many search engines, like MassSpecFile names a format rather than a tool
    }
}
