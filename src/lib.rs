pub mod air;
pub mod error;
mod fiatshamir;
pub mod field;
mod fri;
mod merkle;
pub mod starks;
mod util;

use digest::generic_array::GenericArray;
use digest::OutputSizeUser;

pub type Hash<D> = GenericArray<u8, <D as OutputSizeUser>::OutputSize>;
