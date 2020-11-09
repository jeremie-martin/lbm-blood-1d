#[macro_export]
/// Wrapper for repeated .zip() followed by a .map().collect()
///
/// Example
///
/// ```no_run
/// zip_map!(v.cells, |(A, u, f)| A + u / (f.0 + f.1 + f.2))
/// // =>
/// v.cells.f.iter()
///     .zip(v.cells.u.iter())
///     .zip(v.cells.A.iter())
///     .map(|((f, u), A)| A + u / (f.0 + f.1 + f.2))
/// ```
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

#[macro_export]
/// Wrapper for repeated `.zip()`, adding the prefix `$name` to each `$ident`
///
/// Example
///
/// ```no_run
/// expand_iter!(v.cells, A, u, f)
/// // =>
/// v.cells.f.iter()
///     .zip(v.cells.u.iter())
///     .zip(v.cells.A.iter())
/// ```
macro_rules! expand_iter {
    ( $name: expr, $last:ident ) => {
        $name.$last.iter()
    };

    ( $name: expr, $head:ident $(, $tail:ident)* ) => {
        expand_iter!($name, $($tail),*).zip($name.$head.iter())
    };
}

#[macro_export]
/// Create nested tuples from a list of idents
///
/// Example
///
/// ```no_run
/// expand_names!(A, u, f)
/// // =>
/// ((f, u), A)
/// ```
macro_rules! expand_names {
    ( $last:ident ) => {
        $last
    };

    ( $head:ident $(, $tail:ident)* ) => {
        (expand_names!($($tail),*), $head)
    };
}

#[macro_export]
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

#[macro_export]
macro_rules! expand_iter_mut {
    ( $name: expr, $last:ident ) => {
        (&$name.$last).into_iter()
    };

    ( $name:expr, mut $last:ident ) => {
        (&mut $name.$last).into_iter()
    };

    ( $name: expr, mut $head:ident, $($($tail:ident)+),+) => {
        expand_iter_mut!($name $(, $($tail)+)*).zip(&mut $name.$head)
    };

    ( $name: expr, $head:ident, $($($tail:ident)+),+) => {
        expand_iter_mut!($name $(, $($tail)+)*).zip(&$name.$head)
    };
}

#[macro_export]
macro_rules! expand_names_mut {
    ( $last:ident ) => {
        $last
    };

    ( mut $last:ident ) => {
        $last
    };

    ( mut $head:ident $(, $($tail:ident)+)* ) => {
        (expand_names_mut!($($($tail)+),*), $head)
    };

    ( $head:ident $(, $($tail:ident)+)* ) => {
        (expand_names_mut!($($($tail)+),*), $head)
    };
}

macro_rules! zip_for_each {
    // n = 1
    ( $name:expr, |$head:ident| $body:expr ) => {
        (&mut $name.$head).iter_mut().map(|$head| $body).collect()
    };

    // n > 1
    ( $name:expr, |($($($list:ident)+),+)| $body:expr ) => {
        expand_iter_mut!($name, $($($list)+),+)
        .for_each(|expand_names_mut!($($($list)+),+)| $body )
    };
}
