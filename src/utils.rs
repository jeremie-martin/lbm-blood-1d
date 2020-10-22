macro_rules! expand_names {
    ( $last:ident ) => {
        $last
    };

    ( $head:ident $(, $tail:ident)* ) => {
        (expand_names!($($tail),*), $head)
    };
}

macro_rules! expand_iter {
    ( $name: expr, $last:ident ) => {
        $name.$last.iter()
    };

    ( $name: expr, $head:ident $(, $tail:ident)* ) => {
        expand_iter!($name, $($tail),*).zip($name.$head.iter())
    };
}

macro_rules! zip_map {
    // n = 1
    ( $name:expr, |$head:ident| $body:expr ) => {
        $name.$head.iter().map(|$head| $body ).collect()
    };

    // n > 1
    ( $name:expr, |($head:ident $(, $tail:ident)*)| $body:expr ) => {
        expand_iter!($name, $head $(, $tail)*)
            .map(|expand_names!($head $(, $tail)*)| $body ).collect()
    };
}

macro_rules! zip_map_enumerate {
    // n = 1
    ( $name:expr, |($idx:ident, $head:ident)| $body:expr ) => {
        $name.$head.iter().enumerate().map(|($idx, $head)| $body ).collect()
    };

    // n > 1
    ( $name:expr, |($idx:ident, $head:ident $(, $tail:ident)*)| $body:expr ) => {
        expand_iter!($name, $head $(, $tail)*).enumerate()
            .map(|($idx, expand_names!($head $(, $tail)*))| $body ).collect()
    };
}
