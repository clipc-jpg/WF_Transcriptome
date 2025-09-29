

use std::path::PathBuf;

use dioxus::prelude::*;

use crate::components::*;
use crate::backend::*;


#[component]
pub fn TextInputField(class: String, title: String, data: Signal<String>) -> Element {

    rsx! {
        div {
            class: format!("{} text-input container", &class),
            h3 {
               class: format!("{} text-input title", &class),
               {title}
            }
            input {
                class: format!("{} text-input input-field", &class),
                r#type: "text", // TODO: when run without external launcher, type may be "file"
                value: data,
                oninput: move |event| {
                    data.set(event.data.value().into());
                }
            }
        }
    }
}

#[component]
pub fn NumberInputField(class: String, title: String, data: Signal<i64>) -> Element {

    rsx! {
        div {
            class: format!("{} number-input container", &class),
            h3 {
               class: format!("{} number-input title", &class),
               {title}
            }
            input {
                class: format!("{} number-input input-field", &class),
                r#type: "number", // TODO: when run without external launcher, type may be "file"
                value: data,
                oninput: move |event| {
                    let prev_val: i64 = *data.read();
                    data.set(event.data.value().parse().unwrap_or(prev_val));
                }
            }
        }
    }
}


#[component]
pub fn FileInputField(class: String, title: String, data: SyncSignal<PathBuf>, webserver: Signal<Option<ColonyAuth>>) -> Element {

    rsx! {
        div {
            class: format!("{} file-input container", &class),
            h3 {
               class: format!("{} file-input title", &class),
               {title}
            }
            div {
                class: format!("{} file-input input-field-container", &class),
                input {
                    class: format!("{} file-input input-field", &class),
                    r#type: "text", // TODO: when run without external launcher, type may be "file"
                    value: data().to_string_lossy().to_string(),
                    oninput: move |event| {
                        data.set(event.data.value().into());
                    }
                }
                RequestFilePathButton { class: class.clone(), webserver: webserver, data: data }
            }
        }
    }
}


#[component]
pub fn DirectoryInputField(class: String, title: String, data: SyncSignal<PathBuf>, webserver: Signal<Option<ColonyAuth>>) -> Element {

    let data_string = use_memo(move || data().to_string_lossy().to_string());

    rsx! {
        div {
            class: format!("{} directory-input container", &class),
            h3 {
               class: format!("{} directory-input title", &class),
               {title}
            }
            div {
                class: format!("{} directory-input input-field-container", &class),
                input {
                    class: format!("{} directory-input input-field", &class),
                    r#type: "text", // TODO: when run without external launcher, type may be "file"
                    value: data_string,
                    oninput: move |event| {
                        data.set(event.data.value().into());
                    }
                }
                RequestDirectoryPathButton { class: class.clone(), webserver: webserver, data: data }
            }
        }
    }
}












