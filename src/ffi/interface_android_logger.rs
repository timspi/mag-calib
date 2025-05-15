#[cfg(feature = "java_wrapper")]
use android_logger;
use std::str::FromStr;
// possible log levels: "trace", "debug", info", "warn", "error", "off"
foreign_class!(
    class AndroidLogger {
        fn init_android_logger(level : &str) {
            android_logger::init_once(
                android_logger::Config::default()
                .with_max_level(log::LevelFilter::from_str(level)
                                    .expect("Failed to parse log level"))
                .with_tag("RustMagCalib"),
            );
        }
    }
);