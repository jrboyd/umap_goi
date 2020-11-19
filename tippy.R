tippy_datasets = function(){
    list(
        tippy_this(
            "selDataSource",
            "Choose between different UMAP version."
        ),
        tippy_this(
            "radioItem_SF_original",
            tooltip = "Data as processed by Seth, unclear on details."
        ),
        tippy_this(
            "radioItem_DEGs",
            tooltip = paste("UMAP derived only from genes that are differentially",
                            "expressed in baseline bulk RNA-seq.  Highlights differences",
                            "between wt and df4.")
        ),
        tippy_this(
            "radioItem_refix",
            tooltip = paste("UMAP derived from \"anchor\" genes that are",
                            "in-common between cell types in both wt and df4 backgorunds.",
                            "Highlights cell type differences and allows wt and df4",
                            "cells of to same type to cluster together despite background differences.")
        ),
        tippy_this(
            "radioItem_toy",
            tooltip = "Refixed data filtered for Cd* genes for speed of testing."
        ),
        tippy_this(
            "radioItem_refix_Bcell",
            tooltip = paste("2nd anchored dataset.  Additionally has removed suspected",
                            "ertyhrotcyte, megakaryocyte and mystery island cluster that was df4 biased.")
        )
        # tippy_this(
        #     "radio_SF_original",
        #     tooltip = "Data as processed by Seth, unclear on details."
        # ),
        # tippy_this(
        #     "radio_DEGs",
        #     tooltip = paste("UMAP derived only from genes that are differentially",
        #                     "expressed in baseline bulk RNA-seq.  Highlights differences",
        #                     "between wt and df4.")
        # ),
        # tippy_this(
        #     "radio_refix",
        #     tooltip = paste("UMAP derived from \"anchor\" genes that are",
        #                     "in-common between cell types in both wt and df4 backgorunds.",
        #                     "Highlights cell type differences and allows wt and df4",
        #                     "cells of to same type to cluster together despite background differences.")
        # ),
        # tippy_this(
        #     "radio_toy",
        #     tooltip = "Refixed data filtered for Cd* genes for speed of testing."
        # ),
        # tippy_this(
        #     "radio_refix_Bcell",
        #     tooltip = paste("2nd anchored dataset.  Additionally has removed suspected",
        #                     "ertyhrotcyte, megakaryocyte and mystery island cluster that was df4 biased.")
        # ),
        # tippy_this(
        #     "radioLabel_refix_Bcell",
        #     tooltip = paste("2nd anchored dataset.  Additionally has removed suspected",
        #                     "ertyhrotcyte, megakaryocyte and mystery island cluster that was df4 biased.")
        # ),
        # tippy_this(
        #     "radioLabel_toy",
        #     tooltip = "Refixed data filtered for Cd* genes for speed of testing."
        # ),
        # tippy_this(
        #     "radioLabel_SF_original",
        #     tooltip = "Data as processed by Seth, unclear on details."
        # ),
        # tippy_this(
        #     "radioLabel_DEGs",
        #     tooltip = paste("UMAP derived only from genes that are differentially",
        #                     "expressed in baseline bulk RNA-seq.  Highlights differences",
        #                     "between wt and df4.")
        # ),
        # tippy_this(
        #     "radioLabel_refix",
        #     tooltip = paste("UMAP derived from \"anchor\" genes that are",
        #                     "in-common between cell types in both wt and df4 backgorunds.",
        #                     "Highlights cell type differences and allows wt and df4",
        #                     "cells of to same type to cluster together despite background differences.")
        # ),
        # tippy_this(
        #     "radioLabel_refix_Bcell",
        #     tooltip = paste("2nd anchored dataset.  Additionally has removed suspected",
        #                     "ertyhrotcyte, megakaryocyte and mystery island cluster that was df4 biased.")
        # )
    )
}