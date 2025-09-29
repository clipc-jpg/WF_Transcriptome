

mod components;
mod pages;
mod backend;

use std::io::Read;

use crate::pages::*;
use crate::backend::*;

use dioxus::prelude::*;
use serde::{Deserialize, Serialize};
#[cfg(feature = "web")]
use dioxus_web;

#[cfg(feature = "server")]
use dioxus_fullstack;

#[cfg(feature = "server")]
use dioxus::prelude::DioxusRouterExt;
#[cfg(feature = "server")]
use axum::routing::post;
#[cfg(feature = "server")]
use tokio_stream::{wrappers::WatchStream, StreamExt};

fn main() {

    #[cfg(feature = "web")]
    {
        // Hydrate the application on the client
        //dioxus_web::launch::launch_cfg(App, dioxus_web::Config::new().hydrate(true));
        dioxus_web::launch::launch_cfg(App, dioxus_web::Config::new().hydrate(true));
        //LaunchBuilder::new().launch(App);
    }

     #[cfg(feature = "server")]
     {
         use axum::{response::Html, routing::get, Router};
 		//use dioxus_fullstack::prelude::*;



         tokio::runtime::Runtime::new()
             .unwrap()
             .block_on(async move {

                 unsafe {

                     // termination signal: one route changes a watch channel state (a watch channel sender can be cloned into a route's closure)
                     // this will return indirectly into a tokio::select! statement such that the program will proceed
                     // then, the axum server will be dropped and the program terminates
                     let (transmitter, mut receiver) = tokio::sync::watch::channel(false);

                     // a Watch Stream can be safely dropped inside the tokio::select! statement <=> it implements Unpin
                     //let mut stream = WatchStream::new(receiver).map(|payload| {
 //                        println!("{}", payload);
 //                    });
                     let mut stream = WatchStream::from_changes(receiver);
                     //futures::pin_mut!(stream);

                     // build our application with a route
                     let router = axum::Router::new()
                         // .with_graceful_shutdown(async { receiver.await.ok(); } ) // hidden behind some feature
                         .route("/terminate", get(move || async move { transmitter.send(true); }))
                         // Server side render the application, serve static assets, and register server functions
                         //.serve_dioxus_application("", ServeConfigBuilder::new(App, ()));
                         .serve_dioxus_application(ServeConfigBuilder::new(), App);

                     // run it
                     let listener = tokio::net::TcpListener::bind("127.0.0.1:9283")
                         .await
                         .unwrap();
                     println!("listening on {}", listener.local_addr().unwrap());

                     tokio::select! {
                         _ = axum::serve(listener, router) => {},
                         _ = stream.next() => {println!("Server terminated")}
                     }

                 }

         });
     }
}


//fn main() {
//    #[cfg(feature = "web")]
//    // Hydrate the application on the client
//    dioxus_web::launch::launch_cfg(app, dioxus_web::Config::new().hydrate(true));
//
//    #[cfg(feature = "server")]
//    {
//        use crate::auth::*;
//        use axum::routing::*;
//        use axum_session::SessionConfig;
//        use axum_session::SessionStore;
//        use axum_session_auth::AuthConfig;
//        use axum_session_auth::SessionSqlitePool;
//        simple_logger::SimpleLogger::new().init().unwrap();
//        tokio::runtime::Runtime::new()
//            .unwrap()
//            .block_on(async move {
//                let pool = connect_to_database().await;
//
//                //This Defaults as normal Cookies.
//                //To enable Private cookies for integrity, and authenticity please check the next Example.
//                let session_config = SessionConfig::default().with_table_name("test_table");
//                let auth_config = AuthConfig::<i64>::default().with_anonymous_user_id(Some(1));
//                let session_store = SessionStore::<SessionSqlitePool>::new(
//                    Some(pool.clone().into()),
//                    session_config,
//                )
//                .await
//                .unwrap();
//
//                User::create_user_tables(&pool).await;
//
//                // build our application with some routes
//                let app = Router::new()
//                    // Server side render the application, serve static assets, and register server functions
//                    .serve_dioxus_application(ServeConfig::builder().build(), || {
//                        VirtualDom::new(app)
//                    })
//                    .await
//                    .layer(
//                        axum_session_auth::AuthSessionLayer::<
//                            crate::auth::User,
//                            i64,
//                            axum_session_auth::SessionSqlitePool,
//                            sqlx::SqlitePool,
//                        >::new(Some(pool))
//

































