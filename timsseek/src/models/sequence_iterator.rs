use super::digest::DigestSlice;
use super::query_chunk::NamedQueryChunk;
use crate::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;

pub struct DigestedSequenceIterator {
    digest_sequences: Vec<DigestSlice>,
    chunk_size: usize,
    iteration_index: usize,
    converter: SequenceToElutionGroupConverter,
    build_decoys: bool,
}

impl DigestedSequenceIterator {
    pub fn new(
        digest_sequences: Vec<DigestSlice>,
        chunk_size: usize,
        converter: SequenceToElutionGroupConverter,
        build_decoys: bool,
    ) -> Self {
        Self {
            digest_sequences,
            chunk_size,
            converter,
            iteration_index: 0,
            build_decoys,
        }
    }

    fn get_chunk_digests(&self, chunk_index: usize) -> &[DigestSlice] {
        let start = chunk_index * self.chunk_size;
        if start >= self.digest_sequences.len() {
            return &[];
        }
        let end = start + self.chunk_size;
        let end = if end > self.digest_sequences.len() {
            self.digest_sequences.len()
        } else {
            end
        };
        &self.digest_sequences[start..end]
    }

    fn get_chunk(&self, chunk_index: usize) -> Option<NamedQueryChunk> {
        let seqs = self.get_chunk_digests(chunk_index);
        if seqs.is_empty() {
            return None;
        }
        let (eg_seq, eg_chunk, ei_chunk, charge_chunk) =
            self.converter.convert_sequences(seqs).unwrap();
        let eg_seq = eg_seq.into_iter().cloned().collect();
        Some(NamedQueryChunk::new(
            eg_seq,
            charge_chunk,
            eg_chunk,
            ei_chunk,
        ))
    }

    fn get_decoy_chunk(&self, chunk_index: usize) -> Option<NamedQueryChunk> {
        let seqs = self.get_chunk_digests(chunk_index);
        if seqs.is_empty() {
            return None;
        }
        let decoys = seqs
            .iter()
            .map(|x| x.as_decoy())
            .enumerate()
            .collect::<Vec<(usize, DigestSlice)>>();

        let (eg_seq, eg_chunk, ei_chunk, charge_chunk) = self
            .converter
            .convert_enumerated_sequences(&decoys)
            .unwrap();
        let eg_seq = eg_seq.into_iter().cloned().collect();
        Some(NamedQueryChunk::new(
            eg_seq,
            charge_chunk,
            eg_chunk,
            ei_chunk,
        ))
    }
}

impl Iterator for DigestedSequenceIterator {
    type Item = NamedQueryChunk;

    fn next(&mut self) -> Option<Self::Item> {
        // If its an even iteration, we return the targets.
        // And if its an odd iteration, we return the decoys.
        // IF the struct is requested to build decoys.
        let mut decoy_batch = false;
        let index_use;
        if self.build_decoys {
            index_use = self.iteration_index / 2;
            let decoy_index = self.iteration_index % 2;
            if decoy_index == 1 {
                decoy_batch = true;
            }
            self.iteration_index += 1;
        } else {
            index_use = self.iteration_index;
            self.iteration_index += 1;
        }

        let out = if decoy_batch {
            self.get_decoy_chunk(index_use)
        } else {
            self.get_chunk(index_use)
        };

        if out.is_some() {
            // Make super sure I am not returning an empty chunk
            assert!(!out.as_ref().unwrap().is_empty())
        }

        out
    }
}

impl ExactSizeIterator for DigestedSequenceIterator {
    fn len(&self) -> usize {
        let num_chunks = self.digest_sequences.len() / self.chunk_size;
        if self.build_decoys {
            num_chunks * 2
        } else {
            num_chunks
        }
    }
}
