

use std::net::SocketAddr;
use std::path::PathBuf;
use std::sync::OnceLock;
use std::env;

use dioxus::prelude::*;
use uuid::Uuid;
use serde::{Serialize, Deserialize};
use serde_json;
use regex::Regex;

use crate::pages::*;


static SERVEMODE: OnceLock<Option<ColonyAuth>> = OnceLock::new() as OnceLock<Option<ColonyAuth>>;

#[server]
pub async fn get_serve_mode() -> Result<Option<ColonyAuth>, ServerFnError> {
    // Perform some expensive computation or access a database on the server
    return Ok(SERVEMODE.get_or_init(|| {
        let args: Vec<String> = env::args().into_iter().map(|s| s.trim().to_string()).collect();
        let webserver = if args.contains(&"--colony-interop".to_string()) {
            Some(ColonyAuth::test_example())
        } else {
            None
        };
        webserver
    }).clone());
}


//################################################################################
//## Establishing a connection to Colony Launcher
//################################################################################

#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub enum http_connect {
    http,
    https
}

impl http_connect {
    pub fn as_str(&self) -> &str {
        match self {
            http_connect::https => "https",
            http_connect::http => "http"
        }
    }
}

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct ColonyAuth {
    http_conn: http_connect,
    addr: SocketAddr,
    pw: String
}

impl ColonyAuth {
    //TODO: use Box<dyn Error> as Error content
    pub fn from(s: &str) -> Result<ColonyAuth, ()> {
        let pattern = Regex::new(r#"([^@]+)@(https?\://)(.*)"#).unwrap();
        let (_full,[pw_str, http_string, addr_str]) = pattern.captures(&s).ok_or(())?.extract();

        let http_conn = match http_string {
            "https://" => http_connect::https,
            "http://" => http_connect::http,
            _ => {return Err(());}
        };

        let addr = addr_str.parse().map_err(|_| ())?;

        let pw = match Uuid::parse_str(pw_str) {
            Ok(uuid) => uuid,
            Err(_) => {return Err(());}
        };

        let pw = pw_str.to_string();

        return Ok(ColonyAuth { http_conn, addr, pw });
    }


    pub fn test_example() -> ColonyAuth {
        return ColonyAuth { http_conn: http_connect::http, addr: SocketAddr::from(([127, 0, 0, 1], 20311)) , pw: format!("{:?}", Uuid::new_v4()) };
    }

    pub fn to_address_string(&self) -> String {
        let http_conn = self.http_conn.as_str();
        let addr = self.addr;
        return format!("{http_conn}://{addr}");
    }
}



//################################################################################
//## Transmitting Configuration
//################################################################################

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowConfiguration {
    fastq_directory_path: PathBuf,
    ref_annotation_path: PathBuf,
    ref_genome_path: Option<PathBuf>,
    ref_transcriptome_path: Option<PathBuf>,
    mapping: MappingConfig,
    differential_expression_analysis: Option<DeaConfig>
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MappingConfig {
    indexing: String,
    alignment: String,
    min_mapping_quality: Option<i64>,
    transcriptome_assembly: Option<String>,
}

impl MappingConfig {
    pub fn from_genome(minimap_State: MinimapState) -> MappingConfig {
        return MappingConfig {
            indexing: (minimap_State.indexing)(),
            alignment: (minimap_State.alignment)(),
            min_mapping_quality: Some((minimap_State.min_mapping_quality)()),
            transcriptome_assembly: Some((minimap_State.transcriptome_merge)()),
        };
    }


    pub fn from_trnscriptome(minimap_State: MinimapState_Transcriptome) -> MappingConfig {
        return MappingConfig {
            indexing: (minimap_State.indexing)(),
            alignment: (minimap_State.alignment)(),
            min_mapping_quality: None,
            transcriptome_assembly: None,
        };
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DeaConfig {
    sample_sheet_path: PathBuf,
    minimum_gene_expression: i64,
    minimum_gene_expression_samples: i64,
    minimum_feature_expression: i64,
    minimum_feature_expression_samples: i64,
}

impl DeaConfig {
    pub fn from(dea_state: DeaState) -> DeaConfig {
        return DeaConfig {
            sample_sheet_path: (dea_state.sample_sheet)(),
            minimum_gene_expression: (dea_state.minimum_gene_expression)(),
            minimum_gene_expression_samples: (dea_state.minimum_gene_expression_samples)(),
            minimum_feature_expression: (dea_state.minimum_feature_expression)(),
            minimum_feature_expression_samples: (dea_state.minimum_feature_expression_samples)(),
        };
    }
}


impl WorkflowConfiguration {
    pub fn from(chosen_reference: ChosenReferenceState, dea_state: Option<DeaState>) -> WorkflowConfiguration {


        let differential_expression_analysis = dea_state.map(|state| DeaConfig::from(state));

        match chosen_reference {
            ChosenReferenceState::GenomeReference(genome_ref_state) => {
                let fastq_directory_path = (genome_ref_state.fastq_directory)();
                let ref_annotation_path = (genome_ref_state.ref_annotation)();
                let ref_genome_path = Some((genome_ref_state.ref_genome)());
                let ref_transcriptome_path = None;
                let mapping = MappingConfig::from_genome(genome_ref_state.minimap);

                return WorkflowConfiguration { fastq_directory_path, ref_annotation_path, ref_genome_path, ref_transcriptome_path, mapping, differential_expression_analysis };
            },
            ChosenReferenceState::TranscriptomeReference(transcr_ref_state) => {
                let fastq_directory_path = (transcr_ref_state.fastq_directory)();
                let ref_annotation_path = (transcr_ref_state.ref_annotation)();
                let ref_genome_path = None;
                let ref_transcriptome_path = Some((transcr_ref_state.ref_transcriptome)());
                let mapping = MappingConfig::from_trnscriptome(transcr_ref_state.minimap);

                return WorkflowConfiguration { fastq_directory_path, ref_annotation_path, ref_genome_path, ref_transcriptome_path, mapping, differential_expression_analysis };
            }
        }
    }

    pub fn to_json(&self) -> String {
        return serde_json::to_string_pretty(&self).unwrap_or_else(|_| String::new());
    }
}









