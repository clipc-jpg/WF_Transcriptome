
use std::net::{IpAddr, Ipv4Addr};
use std::path::PathBuf;

use dioxus::prelude::*;


use crate::components::*;
use crate::backend::*;

#[component]
pub fn TestPage() -> Element {
    rsx! {
        "Test Test Test"
    }
}

#[component]
pub fn MainPage(native_webserver: Signal<Option<ColonyAuth>>) -> Element {

    let nucleotide_type = use_signal(move || Nucleotide::cDNA );

    let minimap_indexing = use_signal(move || String::from("-k14") );
    let minimap_alignment = use_signal(move || String::from("-uf") );

    let fastq_dir = use_signal_sync(move || PathBuf::new() );
    let ref_annotation = use_signal_sync(move || PathBuf::new() );
    let genome_reference = ChosenReferenceState::GenomeReference(
        GenomeReferenceState {
            fastq_directory: fastq_dir.clone(),
            ref_annotation: ref_annotation.clone(),
            ref_genome: use_signal_sync(move || PathBuf::new() ),
            minimap: MinimapState {
                indexing: minimap_indexing.clone(),
                alignment: minimap_alignment.clone(),
                min_mapping_quality: use_signal(move || 40 as i64 ),
                transcriptome_merge: use_signal(move || String::from("--conservative") ),
            }
        }
    );

    //let genome_reference = use_signal(move || genome_ref);

    let transcriptome_reference = ChosenReferenceState::TranscriptomeReference(
        TranscriptomeReferenceState {
            fastq_directory: fastq_dir,
            ref_annotation: ref_annotation,
            ref_transcriptome: use_signal_sync(move || PathBuf::new() ),
            minimap: MinimapState_Transcriptome {
                indexing: minimap_indexing,
                alignment: minimap_alignment,
            }
        }
    );

    //let transcriptome_reference = use_signal(move || transcriptome_ref);

    let chosen_reference = use_signal(move || ChosenReference::GenomeReference );
    let chosen_reference_state = use_memo(move || {
        match chosen_reference() {
            ChosenReference::GenomeReference => genome_reference,
            ChosenReference::TranscriptomeReference => transcriptome_reference
        }
    });

    let dea_state = DeaState {
        sample_sheet: use_signal_sync(move || PathBuf::new() ),
        minimum_gene_expression: use_signal(move || 10 as i64 ),
        minimum_gene_expression_samples: use_signal(move || 3 as i64 ),
        minimum_feature_expression: use_signal(move || 3 as i64 ),
        minimum_feature_expression_samples: use_signal(move || 1 as i64 ),
    };

    //let workdir_state = WorkDirState {
//        working_directory: use_signal_sync(move || PathBuf::new() ),
//        output_directory: use_signal_sync(move || PathBuf::new() ),
//    };

    //let dea_state = use_signal(move || dea_state);

    let perform_dea = use_signal(|| true);

    rsx! {
        div {
            class: "main-page responsive-container",
            div {
                class: "main-page background",
                //
                // Header section
                //
                h1 {
                    class: "main-page header title",
                    "WF-Transcriptomes"
                }
                h3 {
                    class: "main-page header subtitle",
                    "an adaptation of "
                    br {}
                    "https://github.com/epi2me-labs/"
                    br {}
                    "wf-transcriptomes"
                }
                //
                // nucleotide section
                //
                div {
                    class: "main-page nucleotide-section section-container",
                    h2 {
                        class: "main-page nucleotide-section title",
                        "Measured Nucleotides"
                    }
                    div {
                        class: "main-page nucleotide-section button-container",
                        ButtonCard  { class: "nucleotide-section cdna-button", selected: nucleotide_type.clone(), choice: Nucleotide::cDNA }
                        ButtonCard  { class: "nucleotide-section rna-button", selected: nucleotide_type.clone(), choice: Nucleotide::RNA }
                    }
                }
                //
                // genome/transcriptome section
                //

                div {
                    class: "main-page reference-section container",
                    div {
                        class: "main-page reference-section choice-tabs",
                        ReferenceChoiceTabs {
                            class: "main-page reference-section",
                            selected: chosen_reference,
                            selected_state: chosen_reference_state,
                            webserver: native_webserver
                        }
                    }
                    div {
                        class: "main-page reference-section choice-content",

                    }
                }

                //
                // differential gene expression analysis section
                //

                DEACard {class: "main-page dea-section", inputs: dea_state, toggle: perform_dea, webserver: native_webserver }

                //
                // final action section
                //

                DownloadConfigJsonButton {
                    class: "main-page".to_string(), webserver: native_webserver,
                    data: (chosen_reference_state, dea_state),
                    perform_dea: perform_dea
                }
            }
        }
    }
}

pub trait ShortDescription {
    fn title(&self) -> &str;
    fn subtitle(&self) -> &str;
    fn description(&self) -> &str;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Nucleotide {
    cDNA,
    RNA
}

impl ShortDescription for Nucleotide {
    fn title(&self) -> &str {
        match self {
            Nucleotide::cDNA => "T",
            Nucleotide::RNA => "U"
        }
    }

    fn subtitle(&self) -> &str {
        match self {
            Nucleotide::cDNA => "cDNA-Experiment",
            Nucleotide::RNA => "RNA-Experiment"
        }
    }

    fn description(&self) -> &str {
        ""
    }
}


#[component]
pub fn ButtonCard<T: 'static + ShortDescription + Clone + PartialEq >(class: String, mut selected: Signal<T>, choice: T) -> Element {
    rsx! {
        div {
            class: if selected() == choice {
                format!("{} button-card button-background selected", &class)
            } else {
                format!("{} button-card button-background", &class)
            },
            div {
                class: if selected() == choice {
                    format!("{} button-card button-container selected", &class)
                } else {
                    format!("{} button-card button-container", &class)
                },
                onclick: move |_| {
                    selected.set(choice.clone());
                },
                h2 {
                    class: format!("{} button-card title", &class),
                    {choice.title()}
                }
                h4 {
                    class: format!("{} button-card subtitle", &class),
                    {choice.subtitle()}
                }
            }
        }
    }
}

pub trait TabPanel {
    fn tab_elements() -> Vec<Element>;
    fn tab_title(&self) -> Element;
    fn subpanel(&self) -> Element;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChosenReference {
    GenomeReference,
    TranscriptomeReference
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChosenReferenceState {
    GenomeReference(GenomeReferenceState),
    TranscriptomeReference(TranscriptomeReferenceState)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenomeReferenceState {
    pub fastq_directory: SyncSignal<PathBuf>,
    pub ref_annotation: SyncSignal<PathBuf>,
    pub ref_genome: SyncSignal<PathBuf>,
    pub minimap: MinimapState,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TranscriptomeReferenceState {
    pub fastq_directory: SyncSignal<PathBuf>,
    pub ref_annotation: SyncSignal<PathBuf>,
    pub ref_transcriptome: SyncSignal<PathBuf>,
    pub minimap: MinimapState_Transcriptome,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MinimapState {
    pub indexing: Signal<String>,
    pub alignment: Signal<String>,
    pub min_mapping_quality: Signal<i64>,
    pub transcriptome_merge: Signal<String>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MinimapState_Transcriptome {
    pub indexing: Signal<String>,
    pub alignment: Signal<String>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeaState {
    pub sample_sheet: SyncSignal<PathBuf>,
    pub minimum_gene_expression: Signal<i64>,
    pub minimum_gene_expression_samples: Signal<i64>,
    pub minimum_feature_expression: Signal<i64>,
    pub minimum_feature_expression_samples: Signal<i64>,
}

//#[derive(Debug, Clone, Copy, PartialEq)]
//pub struct WorkDirState {
//    pub working_directory: SyncSignal<PathBuf>,
//    pub output_directory: SyncSignal<PathBuf>,
//}


impl TabPanel for ChosenReferenceState {
    fn tab_title(&self) -> Element {
        rsx! {
            h2 {
                class: "TabPanel title",
                match &self {
                    &ChosenReferenceState::GenomeReference(_) => "Reference Genome",
                    &ChosenReferenceState::TranscriptomeReference(_) => "Precomputed Transcriptome"
                }
            }
        }
    }

    fn tab_elements() -> Vec<Element> {
        vec! [
            rsx! {
                h2 {
                    class: "TabPanel title",
                    "Reference Genome"
                }
            },
            rsx! {
                h2 {
                    class: "TabPanel title",
                    "Precomputed Transcriptome"
                }
            }
        ]
    }

    fn subpanel(&self) -> Element {
        match &self {
            &ChosenReferenceState::GenomeReference(_) => rsx! {},
            &ChosenReferenceState::TranscriptomeReference(_) => rsx! {}
        }
    }
}

#[component]
pub fn ReferenceChoiceTabs(class: String, mut selected_state: Memo<ChosenReferenceState>, selected: Signal<ChosenReference>, webserver: Signal<Option<ColonyAuth>>) -> Element {
    rsx! {
        div {
            class: format!("{} choice-tab container", &class),
            div {
                class: format!("{} choice-tab header-container", &class),
                div {
                    class: match selected() {
                        ChosenReference::GenomeReference => format!("{} choice-tab genome-reference choice-container selected", &class),
                        ChosenReference::TranscriptomeReference => format!("{} choice-tab genome-reference choice-container", &class),
                    },
                    onclick: move |_| {
                        selected.set(ChosenReference::GenomeReference);
                    },
                    div {
                        class: match selected() {
                            ChosenReference::GenomeReference => format!("{} choice-tab genome-reference choice-area selected", &class),
                            ChosenReference::TranscriptomeReference => format!("{} choice-tab genome-reference choice-area", &class),
                        },
                        h2 {
                            class: format!("{} choice-tab genome-reference element", &class),
                            "Reference Genome"
                        }
                    }
                }
                div {
                    class: match selected() {
                        ChosenReference::GenomeReference => format!("{} choice-tab transcriptome-reference choice-container", &class),
                        ChosenReference::TranscriptomeReference => format!("{} choice-tab transcriptome-reference choice-container selected", &class),
                    },
                    onclick: move |_| {
                        selected.set(ChosenReference::TranscriptomeReference);
                    },
                    div {
                        class: match selected() {
                            ChosenReference::GenomeReference => format!("{} choice-tab transcriptome-reference choice-area", &class),
                            ChosenReference::TranscriptomeReference => format!("{} choice-tab transcriptome-reference choice-area selected", &class),
                        },
                        h2 {
                            class: format!("{} choice-tab transcriptome-reference element", &class),
                            "Precomputed Transcriptome"
                        }
                    }
                }
            }
            match selected_state() {
                ChosenReferenceState::GenomeReference(genome_state) => rsx! { ReferenceGenomeTab { class: class.clone(), inputs: genome_state,webserver } },
                ChosenReferenceState::TranscriptomeReference(transcript_state) => rsx! { PrecomputedTranscriptomeTab {class: class.clone(), inputs: transcript_state, webserver} }
            }
        }
    }
}


#[component]
pub fn ReferenceGenomeTab(class: String, inputs: GenomeReferenceState, webserver: Signal<Option<ColonyAuth>>) -> Element {
    rsx! {
        div {
            class: format!("{} choice-tab genome-reference tab-content", &class),
            h2 {
                class: format!("{} choice-tab genome-reference tab-content-title", &class),
                "General Options"
            }
            DirectoryInputField { class: class.clone(), title: "Fastq Directory", data: inputs.fastq_directory.clone(), webserver }
            FileInputField { class: class.clone(), title: "Reference Genome", data: inputs.ref_genome.clone(), webserver }
            FileInputField { class: class.clone(), title: "Reference Annotation", data: inputs.ref_annotation.clone(), webserver }
            div {
                class: format!("{} choice-tab spacer", &class),
            }
            h2 {
                class: format!("{} choice-tab genome-reference tab-content-title", &class),
                "Mapping Options"
            }
            TextInputField { class: class.clone(), title: "Indexing via minimap2", data: inputs.minimap.indexing.clone() }
            TextInputField { class: class.clone(), title: "Alignment via minimap2", data: inputs.minimap.alignment.clone() }
            NumberInputField { class: class.clone(), title: "Minimum Mapping Quality", data: inputs.minimap.min_mapping_quality.clone() }
            TextInputField { class: class.clone(), title: "Transcriptome Assembly via stringtie2", data: inputs.minimap.transcriptome_merge.clone() }
        }
    }
}

#[component]
pub fn PrecomputedTranscriptomeTab(class: String, inputs: TranscriptomeReferenceState, webserver: Signal<Option<ColonyAuth>>) -> Element{
    rsx! {
        div {
            class: format!("{} choice-tab transcriptome-reference tab-content", &class),
            h2 {
                class: format!("{} choice-tab transcriptome-reference tab-content-title", &class),
                "General Options"
            }
            DirectoryInputField { class: class.clone(), title: "Fastq Directory", data: inputs.fastq_directory.clone(),
            webserver }
            FileInputField { class: class.clone(), title: "Reference Transcriptome", data: inputs.ref_transcriptome.clone(), webserver }
            FileInputField { class: class.clone(), title: "Reference Annotation", data: inputs.ref_annotation.clone(), webserver }
            div {
                class: format!("{} choice-tab spacer", &class)
            }
            h2 {
                class: format!("{} choice-tab genome-reference tab-content-title", &class),
                "Mapping Options"
            }
            TextInputField { class: class.clone(), title: "Indexing via minimap2", data: inputs.minimap.indexing.clone() }
            TextInputField { class: class.clone(), title: "Alignment via minimap2", data: inputs.minimap.alignment.clone() }
        }
    }
}

//################################################################################
//## Input Card for Differential Gene Expression Analysis
//################################################################################

#[component]
pub fn DEACard(class: String, inputs: DeaState, toggle: Signal<bool>, webserver: Signal<Option<ColonyAuth>>) -> Element {
    rsx! {
        div {
            class: format!("{} dea-card container", &class),
            div {
                class: format!("{} dea-card header", &class),
                h2 {
                    class: format!("{} dea-card header-title", &class),
                    "Differential Gene Expression Analysis"
                }
                div {
                    class: format!("{} dea-card toggle-button-container", &class),
                    onclick: move |_| {
                        toggle.toggle();
                    },
                    match toggle() {
                        true => rsx! {
                            svg {
                                width: "16",
                                height: "16",
                                fill: "currentColor",
                                class: "dea-card toggle-button bi bi-x-lg",
                                view_box: "0 0 16 16",
                                path {
                                    d: "M2.146 2.854a.5.5 0 1 1 .708-.708L8 7.293l5.146-5.147a.5.5 0 0 1 .708.708L8.707 8l5.147 5.146a.5.5 0 0 1-.708.708L8 8.707l-5.146 5.147a.5.5 0 0 1-.708-.708L7.293 8z"
                                }
                            }
                        },
                        false => rsx! {
                            svg {
                                width: "16",
                                height: "16",
                                fill: "currentColor",
                                class: "dea-card toggle-button bi bi-x-lg",
                                view_box: "0 0 16 16",
                                path { }
                            }
                        }

                    }
                }
            }
            match toggle() {
                true => rsx! {
                        FileInputField { class: class.clone(), title: "Sample Sheet", data: inputs.sample_sheet.clone(), webserver }
                        NumberInputField { class: class.clone(), title: "Minimum Gene Expression", data: inputs.minimum_gene_expression.clone() }
                        NumberInputField { class: class.clone(), title: "Minimum Samples with Gene Expression", data: inputs.minimum_gene_expression_samples.clone() }
                        NumberInputField { class: class.clone(), title: "Minimum Feature Expression", data: inputs.minimum_feature_expression.clone() }
                        NumberInputField { class: class.clone(), title: "Minimum Samples with Feature Expression", data: inputs.minimum_feature_expression_samples.clone() }
                },
                false => rsx! {div {} }
            }
        }
    }
}



