// allow debuging info in tests
#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn}; // Workaround to use prinltn! for logs.
                                             //
pub mod air;
pub mod fiatshamir;
pub mod field;
pub mod fri;
mod merkle;
pub mod starks;
mod util;

use digest::generic_array::GenericArray;
use digest::OutputSizeUser;

pub type Hash<D> = GenericArray<u8, <D as OutputSizeUser>::OutputSize>;
