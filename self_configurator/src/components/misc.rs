

use dioxus::prelude::*;
use regex::Regex;



pub fn escape(htmlStr: String) -> String {
   return htmlStr
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace("\"", "&quot;")
        .replace("'", "&#39;");

}


pub fn no_escape(htmlStr: String) -> String {
    return htmlStr;
}


pub fn maybe_escape(htmlStr: String) -> String {

    let pattern = Regex::new(r#"[<>"']"#).unwrap();

    if pattern.is_match(&htmlStr) {
        return escape(htmlStr);
    } else {
        return htmlStr;
    }
}




