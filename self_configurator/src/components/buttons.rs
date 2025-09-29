
use std::net::IpAddr;
use std::path::PathBuf;

use dioxus::prelude::*;
//use reqwest::{Client, Response};
use ehttp;

use crate::backend::*;
use crate::ChosenReferenceState;
use crate::DeaState;
//use crate::WorkDirState;
use crate::components::*;


#[component]
pub fn DownloadConfigJsonButton(class: String, webserver: Signal<Option<ColonyAuth>>,
                                data: (Memo<ChosenReferenceState>, DeaState),
                                perform_dea: Signal<bool>) -> Element {

    let config = use_memo( move || {
        let (chosen_ref_state, dea_state) = data;
        let dea_state = if perform_dea() {Some(dea_state)} else {None};
        WorkflowConfiguration::from(chosen_ref_state(), dea_state)
    });

    let config_download = use_memo(move || {
        format!("data:text/plain;charset=utf-8,{}", no_escape(config.read().to_json()))
    });

    rsx! {
        match webserver() {
            Some(ip) => rsx! {
                div {
                    class: "download-button container",
                    //class: class,
                    onclick: move |_| {
                        spawn(async move {
                            let request = ehttp::Request::post("http://127.0.0.1:20311/config/json", config.read().to_json().into_bytes());
                            ehttp::fetch(request, move |result: ehttp::Result<ehttp::Response>| {
                                println!("Response Status code: {:?}", result.clone().map(|resp| resp.status).unwrap_or(404));
                                println!("Response: {:?}", result.clone().map(|resp| {
                                    resp.text().unwrap_or("Empty Rsponse").to_string()
                                }).unwrap_or("Connection Error".to_string()));
                            });
                        });
                    },
                    { "Download configuration" }
                }
            },
            None => rsx! {
                a {
                    class: "download-button",
                    href: config_download,
                    download: "configuration.json",
                    { "Download configuration" }
                }
            }
        }
    }
}


#[component]
pub fn RequestFilePathButton(class: String, webserver: Signal<Option<ColonyAuth>>, data: SyncSignal<PathBuf>) -> Element {

    println!("RequestFilePathButton received {:?}", webserver);
    rsx! {
        div {
            class: match webserver() {
                Some(_) => format!("{} request-file-button container", &class),
                None => format!("{} request-file-button container hidden", &class)
            },
            onclick: move |_| {
                spawn(async move {
                    match webserver() {
                        Some(ip) => {
                            let request = ehttp::Request::get("http://127.0.0.1:20311/choosefile/r1/r2");
                            ehttp::fetch(request, move |result: ehttp::Result<ehttp::Response>| {
                                println!("Response Status code: {:?}", result.clone().map(|resp| resp.status).unwrap_or(404));
                                let resp = result.clone().map(|resp| {
                                    resp.text().unwrap_or("Empty Rsponse").to_string()
                                }).unwrap_or_else(|e| e);
                                println!("Response: {:?}", &resp);
                                data.set(PathBuf::from(resp));
                            });
                        },
                        None => {
                            data.set(PathBuf::from("Missing Connection".to_string()));
                        }
                    }
                });
            },
            svg {
                width: "16",
                height: "16",
                fill: "currentColor",
                class: "request-file-button svg bi bi-file-fill",
                view_box: "0 0 16 16",
                path {
                    fill_rule: "evenodd",
                    d: "M4 0h8a2 2 0 0 1 2 2v12a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2V2a2 2 0 0 1 2-2"
                }
            }
        }
    }
}


#[component]
pub fn RequestDirectoryPathButton(class: String, webserver: Signal<Option<ColonyAuth>>, data: SyncSignal<PathBuf>) -> Element {

    println!("RequestDirectoryPathButton received {:?}", webserver);
    rsx! {
        div {
            class: match webserver() {
                Some(_) => format!("{} request-directory-button container", &class),
                None => format!("{} request-file-button container hidden", &class)
            },
            prevent_default: "onclick",
            prevent_default: "onvalue",
            prevent_default: "submit",
            onclick: move |_| {
                spawn(async move {
                    match webserver() {
                        Some(ip) => {
                            let request = ehttp::Request::get("http://127.0.0.1:20311/choosedirectory/r1/r2");
                            ehttp::fetch(request, move |result: ehttp::Result<ehttp::Response>| {
                                println!("Response Status code: {:?}", result.clone().map(|resp| resp.status).unwrap_or(404));
                                let resp = result.clone().map(|resp| {
                                    resp.text().unwrap_or("Empty Rsponse").to_string()
                                }).unwrap_or_else(|e| e);
                                println!("Response: {:?}", &resp);
                                data.set(PathBuf::from(resp));
                            });
                        },
                        None => {
                            data.set(PathBuf::from("Missing Connection"));
                        }
                    }
                });
            },
            svg {
                width: "16",
                height: "16",
                fill: "currentColor",
                class: match webserver() {
                    Some(_) => format!("{} request-directory-button svg bi bi-file-fill", &class),
                    None => format!("{} request-directory-button svg bi bi-file-fill hidden", &class)
                },
                view_box: "0 0 16 16",
                path {
                    fill_rule: "evenodd",
                    d: "M4 0h8a2 2 0 0 1 2 2v12a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2V2a2 2 0 0 1 2-2"
                }
            }
        }
    }
}









