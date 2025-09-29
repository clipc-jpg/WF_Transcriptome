
use std::env;

use dioxus::prelude::*;

use crate::pages::*;
use crate::components::*;
use crate::backend::*;



#[component]
pub fn App() -> Element {

    let webserver = use_server_future(get_serve_mode)?.value();
    //let webserver = use_signal(|| None as Option<Result<Option<ColonyAuth>, ()>>);

    let mut webserver = use_signal(|| {
        match webserver() {
            Some(Ok(serv_opt)) => serv_opt,
            _ => None
        }
    });

    rsx! {
        style { { include_str!("./../public/norm.css") } }
        style { { include_str!("./../public/resets.css") } }
        style { { include_str!("./../public/styles.css") } }
        //style { { include_str!("./../public/color_palette_green1.css") } }
        style { { include_str!("./../public/color_palette_green2.css") } } // vars replaced with values by hand
        MainPage { native_webserver: webserver }
    }
}


















