"""Highlighter plot tools for sequence alignment visualization.

This module provides classes and utilities for identifying and visualizing
mutations, mismatches, and sequence features in multiple sequence alignments.
It extends the ``Bio.Align.AlignInfo`` and ``Bio.Graphics`` modules with
highlighter-style plots commonly used in HIV/viral sequence analysis.

Classes
-------
Highlighter
    Extracts mismatch and match information from a multiple sequence alignment.
HighlighterPlot
    Renders alignment comparison plots as SVG or PDF using ReportLab.

Functions
---------
codon_position
    Determines the codon position of a residue, accounting for alignment gaps.

Examples
--------
Basic mismatch analysis::

    from Bio import AlignIO
    from Bio.Align import AlignInfo

    alignment = AlignIO.read("sequences.fasta", "fasta")
    highlighter = AlignInfo.Highlighter(alignment, seq_type="NT")
    mismatches = highlighter.list_mismatches(references=0, apobec=True)

Generating a highlighter plot::

    from Bio import Graphics

    plot = Graphics.HighlighterPlot(alignment, seq_type="NT")
    plot.draw_mismatches("output.svg", reference=0, apobec=True)

Notes
-----
This module monkey-patches ``Bio.Align.AlignInfo``, ``Bio.Graphics``, and
``Bio.SeqUtils`` to register its classes and functions into the Biopython
namespace (e.g. ``AlignInfo.Highlighter``, ``Graphics.HighlighterPlot``,
``SeqUtils.codon_position``).
"""

from functools import cache
from typing import Union

import Bio
from Bio import Graphics
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Highlighter:
    """Extract mismatch and match information from a multiple sequence alignment.

    This class provides methods to identify positions in an alignment that
    differ from a reference sequence, and to classify those differences by
    type (e.g. nucleotide substitution, APOBEC signature, stop codon, or
    N-linked glycosylation site). Results can be exported as structured text
    or passed directly to :class:`HighlighterPlot` for visualization.

    Parameters
    ----------
    alignment : list-like of SeqRecord
        A multiple sequence alignment, typically produced by
        ``Bio.AlignIO.read``.
    seq_type : {'NT', 'AA'}
        Sequence type. ``'NT'`` for nucleotide, ``'AA'`` for amino acid.
        This parameter is required.
    codon_offset : int, optional
        Offset to apply when computing codon positions (default 0). Only
        relevant for nucleotide alignments. The value is taken modulo 3.

    Raises
    ------
    ValueError
        If ``seq_type`` is not ``'NT'`` or ``'AA'``.

    Examples
    --------
    >>> from Bio import AlignIO
    >>> from Bio.Align import AlignInfo
    >>> alignment = AlignIO.read("sequences.fasta", "fasta")
    >>> h = AlignInfo.Highlighter(alignment, seq_type="NT")
    >>> mismatches = h.list_mismatches(references=0)
    """

    def __init__(self, alignment, *, seq_type: str=None, codon_offset: int=0):
        """Initialize the Highlighter object.

        Parameters
        ----------
        alignment : list-like of SeqRecord
            A multiple sequence alignment.
        seq_type : {'NT', 'AA'}
            Sequence type. ``'NT'`` for nucleotide, ``'AA'`` for amino acid.
        codon_offset : int, optional
            Reading frame offset (default 0). Taken modulo 3.
        """

        self.alignment = alignment
        self.codon_offet: int=codon_offset % 3

        if seq_type not in ("NT", "AA"):
            raise ValueError(f"type must be provided (either 'NT' or 'AA', got: '{seq_type}')")
        else:
            self.seq_type = seq_type
        
    def get_seq_index_by_id(self, id: str) -> int:
        """Return the index of a sequence in the alignment by its identifier.

        Parameters
        ----------
        id : str
            The sequence identifier to search for (matches ``SeqRecord.id``).

        Returns
        -------
        int
            Zero-based index of the matching sequence in the alignment.

        Raises
        ------
        IndexError
            If no sequence with the given identifier is found.
        """

        for index, sequence in enumerate(self.alignment):
            if sequence.id == id:
                return index
        
        raise IndexError(f"Could not find sequence with id {id}")

    def list_mismatches(self, *, references: Union[int, str]=0, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, codon_offset: int=0) -> list[dict[str: list]]:
        """Return per-sequence mismatch dictionaries relative to a reference.

        Iterates over every sequence in the alignment and identifies positions
        that differ from the designated reference sequence. Optionally
        annotates biologically relevant mutation types.

        Parameters
        ----------
        references : int or str, optional
            The reference sequence, specified either as a zero-based alignment
            index (int) or a sequence identifier (str). Defaults to ``0``
            (the first sequence).
        apobec : bool, optional
            If ``True``, flag G→A substitutions consistent with APOBEC
            editing (GA context, not followed by C). Only applies to
            nucleotide alignments (default ``False``).
        g_to_a : bool, optional
            If ``True``, flag all G→A substitutions regardless of context.
            Only applies to nucleotide alignments (default ``False``).
        stop_codons : bool, optional
            If ``True``, identify positions that introduce a stop codon
            (TAA, TAG, or TGA). Only applies to nucleotide alignments
            (default ``False``).
        glycosylation : bool, optional
            If ``True``, identify asparagine (N) residues in N-X-S/T sequons
            that constitute potential N-linked glycosylation sites. Only
            applies to amino acid alignments (default ``False``).
        codon_offset : int, optional
            Reading frame offset applied when detecting stop codons
            (default ``0``).

        Returns
        -------
        list of dict
            One dictionary per sequence in the alignment. Each dictionary
            maps zero-based alignment positions (int) to a list of
            annotation strings. The first element of each list is the
            residue at that position (or ``'Gap'`` for ``'-'``); subsequent
            elements are one or more of ``'G->A mutation'``, ``'APOBEC'``,
            ``'Stop codon'``, or ``'Glycosylation'``.

        Raises
        ------
        TypeError
            If ``references`` is not an int or str.
        ValueError
            If the reference sequence is empty.
        """

        mismatches: list[dict[str: list]] = []

        if isinstance(references, str):
            references = self.get_seq_index_by_id(references)
        elif not isinstance(references, int):
            raise TypeError(f"Expected reference to be an int or str, got {type(references)}")

        self.references = references

        reference_object = self.alignment[references]

        if not reference_object:
            raise ValueError(f"Reference sequence {references} is empty")
        
        for sequence in self.alignment:
            mismatches.append(self.get_mismatches(sequence=sequence, references=reference_object, seq_type=self.seq_type, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=codon_offset))

        return mismatches

    @staticmethod
    def get_mismatches(*, sequence: Union[str, Seq, SeqRecord], references: Union[str, Seq, SeqRecord], seq_type: str=None, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, codon_offset: int=0) -> dict[int: list]:
        """Return a mismatch dictionary comparing one sequence to a reference.

        Accepts ``str``, :class:`~Bio.Seq.Seq`, or
        :class:`~Bio.SeqRecord.SeqRecord` inputs and normalises them before
        delegating to :meth:`get_mismatches_from_str` (which is cached).

        Parameters
        ----------
        sequence : str, Seq, or SeqRecord
            The query sequence to compare against the reference.
        references : str, Seq, or SeqRecord
            The reference sequence. Must be the same length as ``sequence``.
        seq_type : {'NT', 'AA'}
            Sequence type. Required.
        apobec : bool, optional
            Flag APOBEC-signature G→A mutations (default ``False``).
        g_to_a : bool, optional
            Flag all G→A mutations (default ``False``).
        stop_codons : bool, optional
            Flag positions that introduce stop codons (default ``False``).
        glycosylation : bool, optional
            Flag N-linked glycosylation sequons (default ``False``).
        codon_offset : int, optional
            Reading frame offset for stop codon detection (default ``0``).

        Returns
        -------
        dict
            Maps zero-based alignment positions to annotation lists.
            See :meth:`list_mismatches` for the list format.

        Raises
        ------
        TypeError
            If ``sequence`` or ``references`` is not a str, Seq, or SeqRecord.
        ValueError
            If ``seq_type`` is invalid, or if the sequences differ in length.
        """

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")

        if isinstance(sequence, Seq):
                sequence = str(sequence)
        elif isinstance(sequence, SeqRecord):
                sequence = str(sequence.seq)
        elif not isinstance(sequence, str):
                raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")

        if isinstance(references, Seq):
                references = str(references)
        elif isinstance(references, SeqRecord):
                references = str(references.seq)
        elif not isinstance(references, str):
                raise TypeError(f"Expected reference to be a string, Seq, or SeqRecord, got {type(references)}")
            
        if len(sequence) != len(references):
            raise ValueError("Reference and sequence must be the same length")
        
        return Highlighter.get_mismatches_from_str(sequence=sequence, references=references, seq_type=seq_type, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=codon_offset)
    
    @staticmethod
    @cache
    def get_mismatches_from_str(*, sequence: str, references: str, seq_type: str, apobec: bool, g_to_a: bool, stop_codons: bool=False, glycosylation: bool, codon_offset: int=0) -> dict[int: list]:
        """Return a mismatch dictionary from plain string inputs (cached).

        This method contains the core comparison logic and is separated from
        :meth:`get_mismatches` so that results can be cached using
        :func:`functools.cache`. :class:`~Bio.Seq.Seq` and
        :class:`~Bio.SeqRecord.SeqRecord` objects are not hashable and
        therefore cannot be used as cache keys directly.

        Parameters
        ----------
        sequence : str
            The query sequence (gaps represented as ``'-'``).
        references : str
            The reference sequence. Must be the same length as ``sequence``.
        seq_type : {'NT', 'AA'}
            Sequence type.
        apobec : bool
            Flag APOBEC-signature G→A mutations.
        g_to_a : bool
            Flag all G→A mutations.
        stop_codons : bool, optional
            Flag positions that introduce stop codons (default ``False``).
        glycosylation : bool
            Flag N-linked glycosylation sequons.
        codon_offset : int, optional
            Reading frame offset for stop codon detection (default ``0``).

        Returns
        -------
        dict
            Maps zero-based alignment positions (int) to annotation lists.
            Positions identical to the reference are omitted. Each list
            begins with the residue character (or ``'Gap'``) and may include
            any of ``'G->A mutation'``, ``'APOBEC'``, ``'Stop codon'``, or
            ``'Glycosylation'``.

        Notes
        -----
        Newline characters are stripped from both sequences before comparison.
        Stop codon detection uses :func:`Bio.SeqUtils.codon_position` to
        handle gap-adjusted reading frames. Glycosylation detection follows
        the N-X-S/T sequon rule (where X is not proline).
        """

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")

        sequence = sequence.replace("\n", "")
        references = references.replace("\n", "")

        mismatches: dict = {}

        if sequence == references and (seq_type == "NT" and not stop_codons) and (seq_type == "AA" and not glycosylation):
            return mismatches

        for base_index in range(len(sequence)):
            if references[base_index] != sequence[base_index]:
                mismatches[base_index] = []
                
                if sequence[base_index] != "-":
                    mismatches[base_index].append(sequence[base_index])
                else:
                    mismatches[base_index].append("Gap")

                # APOBEC and G->A mutations only apply to NT sequences
                if seq_type == "NT":
                    if references[base_index] == "G" and sequence[base_index] == "A":
                        if g_to_a:
                            mismatches[base_index].append("G->A mutation")

                        if apobec and base_index <= len(sequence)-3 and sequence[base_index+1] in "AG" and sequence[base_index+2] != "C":
                            mismatches[base_index].append("APOBEC")

            # Stop codons only apply to NT sequences
            if seq_type == "NT":
                if stop_codons and sequence[base_index] in "TU" and base_index <= len(sequence)-3 and SeqUtils.codon_position(sequence, base_index, codon_offset=codon_offset) == 0:
                    base_snippet: str = ""
                    snippet_index: int = base_index+1

                    while len(base_snippet) < 2 and snippet_index < len(sequence):
                        base_snippet += sequence[snippet_index] if sequence[snippet_index] != "-" else ""
                        snippet_index += 1

                        if base_snippet in ("AA", "AG", "GA"):
                            if base_index not in mismatches:
                                mismatches[base_index] = []

                            mismatches[base_index].append("Stop codon")

            # Glycosylation only applies to AA sequences
            elif seq_type == "AA":
                if glycosylation and sequence[base_index] == "N" and base_index <= len(sequence)-3:
                    base_snippet: str = ""
                    snippet_index: int = base_index+1

                    while len(base_snippet) < 2 and snippet_index < len(sequence):
                        base_snippet += sequence[snippet_index] if sequence[snippet_index] != "-" else ""
                        snippet_index += 1
                    
                    if base_snippet[0] != "P" and base_snippet[1] in "ST":
                        if base_index not in mismatches:
                            mismatches[base_index] = []

                        mismatches[base_index].append("Glycosylation")
        
        return mismatches
    
    def export_mismatches(self, output_file, *, references: Union[int, str]=0, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, codon_offset: int=0) -> None:
        """Write mismatch annotations for all sequences to a text file.

        For each sequence in the alignment the output lists the sequence
        identifier followed by one line per annotation type giving the
        (1-based) alignment positions where that type was observed.

        Example output::

            seq1
            A [12 45 67]
            Gap [103]

            seq2
            C [8]

        Parameters
        ----------
        output_file : str or path-like
            Destination file path. The file is written in text mode and
            will be created or overwritten.
        references : int or str, optional
            Reference sequence index or identifier (default ``0``).
        apobec : bool, optional
            Flag APOBEC-signature G→A mutations (default ``False``).
        g_to_a : bool, optional
            Flag all G→A mutations (default ``False``).
        stop_codons : bool, optional
            Flag stop-codon-introducing positions (default ``False``).
        glycosylation : bool, optional
            Flag N-linked glycosylation sequons (default ``False``).
        codon_offset : int, optional
            Reading frame offset (default ``0``).
        """

        output: str = ""
        mismatches = self.list_mismatches(references=references, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=codon_offset)

        for sequence_index, sequence in enumerate(self.alignment):
            working: dict = {}

            for base, codes in mismatches[sequence_index].items():
                for code in codes:
                    if code not in working:
                        working[code] = []
                    
                    working[code].append(base+1)
            
            output += f"{sequence.id}\n"

            for code in sorted(working, key=lambda x: (len(x), x)):
                output += f"{code} ["
                output += " ".join([str(thing) for thing in working[code]])
                output += "]\n"

            output += "\n"

        with open(output_file, mode="wt") as file:
            file.write(output)

    def list_matches(self, *, references=0) -> list[dict[str: list]]:
        """Return per-sequence match dictionaries relative to one or more references.

        For each sequence, identifies positions that *match* at least one
        reference and classifies them by which reference(s) they match.
        Positions that match all references simultaneously are omitted
        (they carry no discriminating information).

        Parameters
        ----------
        references : int, str, Seq, or list thereof, optional
            One or more reference sequences. Each entry may be:

            - **int** – zero-based index into the alignment.
            - **str** – sequence identifier (matched against ``SeqRecord.id``).
            - **Seq** – a raw sequence object.

            A single non-list value is wrapped in a list automatically.
            Defaults to ``0`` (the first sequence).

        Returns
        -------
        list of dict
            One dictionary per sequence in the alignment. For reference
            sequences the dictionary is empty. For all other sequences,
            keys are zero-based alignment positions and values are lists
            containing the zero-based index of every reference that shares
            the same residue at that position, or ``['Unique']`` if no
            reference matches.

        Raises
        ------
        IndexError
            If an integer reference index is out of range.
        TypeError
            If a reference entry is not an int, str, Seq, or dict.
        """

        matches: list[dict[str: list]] = []

        if not isinstance(references, list):
            references = [references]

        reference_objects: list = []
        self.references = []

        for reference in references:
            if isinstance(reference, str):
                self.references.append(self.get_seq_index_by_id(reference).replace("\n", ""))

            elif isinstance(reference, int):
                if reference >= len(self.alignment):
                    raise IndexError(f"Reference index {reference} is out of range")
                
                self.references.append(reference)

            elif isinstance(reference, Seq):
                self.references.append(str(reference).replace("\n", ""))

            elif isinstance(reference, dict):
                if "sequence" not in reference or "color" not in reference:
                    raise ValueError("Reference dictionaries must contain 'sequence' and 'color' keys")

            else:
                raise TypeError(f"Expected reference to be an int or str, got {type(reference)}")
            
            if isinstance(reference, Seq):
                reference_objects.append(str(reference))
            else:
                reference_objects.append(self.alignment[self.references[-1]])

        for sequence_index, sequence in enumerate(self.alignment):
            if sequence_index in self.references:
                matches.append({})
            else:
                matches.append(self.get_matches(sequence=sequence, references=reference_objects, seq_type=self.seq_type))

        return matches
    
    @staticmethod
    def get_matches(*, sequence: Union[str, Seq, SeqRecord], references: Union[list[str, Seq, SeqRecord], str, Seq, SeqRecord], seq_type: str=None) -> dict[int: list]:
        """Return a match dictionary comparing one sequence to one or more references.

        Accepts ``str``, :class:`~Bio.Seq.Seq`, or
        :class:`~Bio.SeqRecord.SeqRecord` inputs, normalises them to plain
        strings, then delegates to :meth:`get_matches_from_str`.

        Parameters
        ----------
        sequence : str, Seq, or SeqRecord
            The query sequence.
        references : str, Seq, SeqRecord, or list thereof
            One or more reference sequences. All must be the same length as
            ``sequence``. A single non-list value is wrapped automatically.
        seq_type : {'NT', 'AA'}
            Sequence type. Required.

        Returns
        -------
        dict
            Maps zero-based alignment positions to match annotation lists.
            See :meth:`list_matches` for the list format. Positions where
            the query matches *all* references are omitted.

        Raises
        ------
        TypeError
            If ``sequence`` or any reference is not a str, Seq, or SeqRecord.
        ValueError
            If ``seq_type`` is invalid or sequences differ in length.
        """
        
        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")
        
        if isinstance(sequence, Seq):
                sequence = str(sequence).replace("\n", "")
        elif isinstance(sequence, SeqRecord):
                sequence = str(sequence.seq).replace("\n", "")
        elif not isinstance(sequence, str):
                raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")
        
        if not isinstance(references, list):
            references = [references]

        new_references: list = []

        for reference in references:
            if isinstance(reference, Seq):
                reference = str(reference).replace("\n", "")
            elif isinstance(reference, SeqRecord):
                reference = str(reference.seq).replace("\n", "")
            elif isinstance(reference, str):
                reference = reference.replace("\n", "")
            elif not isinstance(reference, str):
                raise TypeError(f"Expected reference to be a string, Seq, or SeqRecord, got {type(reference)}")
            
            if len(sequence) != len(reference):
                raise ValueError("All references and sequence must be the same length")
            
            new_references.append(reference)

        return Highlighter.get_matches_from_str(sequence=sequence, references=tuple(new_references), seq_type=seq_type)

    @staticmethod
    @cache
    def get_matches_from_str(*, sequence: str, references: tuple[str], seq_type: str) -> dict[int: list]:
        """Return a match dictionary from plain string inputs (cached).

        Separated from :meth:`get_matches` so results can be cached via
        :func:`functools.cache`. References are passed as a tuple (rather
        than a list) because tuples are hashable.

        Parameters
        ----------
        sequence : str
            The query sequence.
        references : tuple of str
            One or more reference sequences, all the same length as
            ``sequence``.
        seq_type : {'NT', 'AA'}
            Sequence type.

        Returns
        -------
        dict
            Maps zero-based alignment positions to lists. Each list contains
            the zero-based indices of references that match the query at that
            position, or ``['Unique']`` if no reference matches. Positions
            where the query matches *all* references are omitted. The ``'X'``
            residue in a reference is treated as a wildcard.
        """

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")

        matches: dict = {}
        
        # if sequence in references:
        #     return matches

        for base_index in range(len(sequence)):
            matches[base_index] = []
            for reference_index, reference in enumerate(references):
                if reference[base_index] == sequence[base_index] or reference[base_index] == "X":
                    matches[base_index].append(reference_index)

            if not matches[base_index]:
                matches[base_index].append("Unique")

            elif len(matches[base_index]) == len(references):
                del matches[base_index]

        return matches
    
    def export_matches(self, output_file, *, references: tuple[Union[int, str]]=0) -> None:
        """Write match annotations for all sequences to a text file.

        For each sequence the output lists the sequence identifier (with a
        reference label where applicable) followed by one line per match
        category giving the 1-based alignment positions.

        Match categories:

        - **R1, R2, …** – positions matching only that reference.
        - **Multiple matches** – positions matching more than one reference.
        - **Unique in query** – positions matching no reference.

        Example output::

            ref_seq (R1)
            R1 [1 2 3 ...]

            query_seq
            R1 [5 12 20]
            Unique in query [8 15]

        Parameters
        ----------
        output_file : str or path-like
            Destination file path. Written in text mode; created or
            overwritten.
        references : int, str, or tuple thereof, optional
            Reference sequence(s) by index or identifier (default ``0``).
        """

        output: str = ""
        matches = self.list_matches(references=references)
        
        for sequence_index, sequence in enumerate(self.alignment):
            matched: dict = {}

            for base, codes in matches[sequence_index].items():
                if "Unique" in codes:
                    if "Unique" not in matched:
                        matched["Unique"] = []
                    
                    matched["Unique"].append(base+1)
                
                elif len(codes) > 1:
                    if "Multiple" not in matched:
                        matched["Multiple"] = []
                    
                    matched["Multiple"].append(base+1)

                else:
                    if codes[0] not in matched:
                        matched[codes[0]] = []
                    
                    matched[codes[0]].append(base+1)

            if sequence_index in self.references:
                output += f"{sequence.id} (R{self.references.index(sequence_index)+1})\n"
            else:
                output += f"{sequence.id}\n"

            for code in list(range(10)) + ["Unique", "Multiple"]:
                if code not in matched:
                    continue

                if code == "Unique":
                    output += f"Unique in query "

                elif code == "Multiple":
                    output += f"Multiple matches "
                    
                else:
                    output += f"R{code+1} "

                output += f"[{' '.join([str(thing) for thing in matched[code]])}]\n"

            output += "\n"

        with open(output_file, mode="wt") as file:
            file.write(output)
    
AlignInfo.Highlighter = Highlighter


from reportlab.lib import colors
from reportlab.lib.colors import Color
from reportlab.lib.units import inch

from reportlab.pdfbase.pdfmetrics import stringWidth

from reportlab.graphics.shapes import Drawing, String, Line, Rect, Circle, PolyLine, Polygon

from Bio import SeqUtils
from Bio.Graphics import _write
from Bio.Align import AlignInfo


class HighlighterPlot:
    """Render alignment comparison plots as SVG or PDF.

    Produces highlighter-style figures in which each sequence is drawn as a
    horizontal line and coloured marks indicate positions that differ from
    (or match) a designated reference sequence. The visual style is modelled
    on the Los Alamos National Laboratory (LANL) HIV sequence database
    highlighter tool.

    Two plot modes are supported:

    - **Mismatch plot** – marks positions where a sequence differs from a
      single reference, coloured by nucleotide or amino acid identity.
    - **Match plot** – marks positions where a sequence matches one or more
      references, coloured by which reference(s) are matched.

    Parameters
    ----------
    alignment : list-like of SeqRecord
        A multiple sequence alignment.
    seq_type : {'NT', 'AA'}, optional
        Sequence type. If not provided (or not ``'NT'``/``'AA'``), the type
        is inferred automatically via :meth:`guess_alignment_type`.
    tree : Bio.Phylo.BaseTree.Tree, optional
        A phylogenetic tree used to order sequences when ``sort='tree'`` is
        passed to a draw method.
    plot_width : float, optional
        Width of the sequence plot area in ReportLab units
        (default ``4 * inch``).
    seq_name_font : str, optional
        Font name for sequence labels (default ``'Helvetica'``).
    seq_name_font_size : int, optional
        Font size for sequence labels in points (default ``8``).
    seq_gap : int or None, optional
        Vertical gap between sequence tracks in ReportLab units. If
        ``None``, defaults to one-fifth of the track height.
    left_margin : float, optional
        Left margin in ReportLab units (default ``0.25 * inch``).
    top_margin : float, optional
        Top margin in ReportLab units (default ``0.25 * inch``).
    bottom_margin : float, optional
        Bottom margin in ReportLab units (default ``0``).
    right_margin : float, optional
        Right margin in ReportLab units (default ``0``).
    plot_label_gap : float, optional
        Gap between the right edge of the plot and the sequence labels
        (default ``inch / 4``).
    mark_reference : bool, optional
        If ``True``, append a reference indicator (e.g. ``(r)``) to the
        label of the reference sequence (default ``True``).
    title_font : str, optional
        Font name for the plot title (default ``'Helvetica'``).
    title_font_size : int, optional
        Font size for the plot title in points (default ``12``).
    ruler : bool, optional
        If ``True``, draw a ruler along the bottom of the plot
        (default ``True``).
    ruler_font : str, optional
        Font name for ruler labels (default ``'Helvetica'``).
    ruler_font_size : int, optional
        Font size for ruler labels in points (default ``6``).
    ruler_major_ticks : int, optional
        Number of major tick marks on the ruler (default ``10``).
    ruler_minor_ticks : int, optional
        Number of minor ticks between each pair of major ticks (default ``3``).
    codon_offset : int, optional
        Reading frame offset for codon-based annotations (default ``0``).

    Raises
    ------
    TypeError
        If ``tree`` is not a ``Bio.Phylo.BaseTree.Tree`` instance.

    Examples
    --------
    >>> from Bio import AlignIO, Graphics
    >>> from reportlab.lib.units import inch
    >>> alignment = AlignIO.read("sequences.fasta", "fasta")
    >>> plot = Graphics.HighlighterPlot(alignment, seq_type="NT", plot_width=6*inch)
    >>> plot.draw_mismatches("out.svg", reference=0, apobec=True)
    """

    def __init__(self, alignment, *, seq_type: str=None, tree: Union[str, object]=None, plot_width: int = 4*inch, seq_name_font: str="Helvetica", seq_name_font_size: int=8, seq_gap: int=None, left_margin: float=.25*inch, top_margin: float=.25*inch, bottom_margin: float=0, right_margin: float=0, plot_label_gap: float=(inch/4), mark_reference: bool=True, title_font="Helvetica", title_font_size: int=12, ruler: bool=True, ruler_font: str="Helvetica", ruler_font_size: int=6, ruler_major_ticks: int=10, ruler_minor_ticks=3, codon_offset: int=0):
        """Initialize the HighlighterPlot object.

        See class docstring for parameter descriptions.
        """

        self.alignment = alignment
        self.codon_offset = codon_offset % 3

        if seq_type not in ("NT", "AA"):
            self.seq_type = self.guess_alignment_type(alignment)
        else:
            self.seq_type = seq_type

        if tree is not None:
            if isinstance(tree, Bio.Phylo.BaseTree.Tree):
                self.tree = tree
            else:
                raise TypeError("Tree must be a Bio.Phylo.BaseTree.Tree object (or a derivative)")

        self._seq_count = len(alignment)
        self._seq_length = len(alignment[0])
        self.mark_reference: bool = mark_reference

        self.seq_name_font: str = seq_name_font
        self.seq_name_font_size: int = seq_name_font_size

        self.top_margin: float = top_margin
        self.bottom_margin: float = bottom_margin
        self.left_margin: float = left_margin
        self.right_margin: float = right_margin
        self.plot_label_gap: float = plot_label_gap

        self.title_font: str = title_font
        self.title_font_size: int = title_font_size

        self.ruler: bool = ruler
        self.ruler_font: str = ruler_font
        self.ruler_font_size: int = ruler_font_size
        self.ruler_major_ticks: int = ruler_major_ticks
        self.ruler_minor_ticks: int = ruler_minor_ticks
        self._ruler_font_height: float = self.ruler_font_size

        self.plot_width: float = plot_width

        self._seq_height: float = self.seq_name_font_size
        self.seq_gap: float = self._seq_height / 5 if seq_gap is None else seq_gap

        self.match_plot_colors: dict[str: dict] = {
            "ML": {
                "references": ["#FF0000", "#537EFF", "#00CB85", "#000000", "#FFA500"],
                "unique": "#EFE645",
                "multiple": "#808080",
            },
            "LANL": {
                "references": ["#ED1C24", "#235192", "#FFC20E", "#00A651", "#8DC73F"],
                "unique": "#000000",
                "multiple": "#666666",
            }
        }

        self.mismatch_plot_colors: dict[str: [dict[str: str]]] = {
            "NT": {
                "LANL": {
                    "A": "#42FF00", 
                    "C": "#41B8EE", 
                    "G": "#FFA500", 
                    "T": "#EE0B10", 
                    "Gap": "#666666",
                },
                "ML": {
                    "A": "#36b809", 
                    "C": "#1282b5", 
                    "G": "#FFA500", 
                    "T": "#c48002", 
                    "Gap": "#666666",
                }
            },
            "AA": {
                "LANL": {
                    "H": "#FF0000",
                    "D": "#302ECD",
                    "E": "#302ECD",
                    "K": "#23659B",
                    "N": "#23659B",
                    "Q": "#23659B",
                    "R": "#23659B",
                    "M": "#2F9A2F",
                    "I": "#42FF00",
                    "L": "#42FF00",
                    "V": "#42FF00",
                    "F": "#F900FF",
                    "W": "#F900FF",
                    "Y": "#F900FF",
                    "C": "#CD2F2E",
                    "A": "#F9CE2E", 
                    "G": "#F9CE2E",
                    "S": "#F9CE2E",
                    "T": "#F9CE2E",
                    "P": "#FBFF00",
                    "Other": "#000000",
                    "Gap": "#bebebe",
                }
            }
        }

    def _setup_drawing(self, *, plot_type: str, output_format: str="svg", title: str=None, sort: str="similar", mark_width: float=1, scale: float=1, sequence_labels: bool=True, legend_entries: list=None):
        """Configure and initialise the ReportLab Drawing object.

        Calculates overall figure dimensions, creates the
        :class:`~reportlab.graphics.shapes.Drawing`, draws the title,
        legend (if entries are provided), and ruler (if enabled), adds
        sequence baselines, and resolves the display order of sequences.

        This method is called internally by :meth:`draw_mismatches` and
        :meth:`draw_matches` before marks are rendered.

        Parameters
        ----------
        plot_type : {'mismatch', 'match'}
            Controls baseline colouring and reference-label logic.
        output_format : str, optional
            Output file format passed to ``Bio.Graphics._write``
            (default ``'svg'``).
        title : str or None, optional
            Plot title drawn above the figure (default ``None``).
        sort : {'similar', 'tree', 'none'}, optional
            Sequence ordering strategy.  ``'similar'`` groups sequences by
            similarity to the reference; ``'tree'`` uses the order of
            terminals in the supplied phylogenetic tree; any other value
            preserves the original alignment order (default ``'similar'``).
        mark_width : float, optional
            Fractional width of each mark relative to one alignment column
            (default ``1``).
        scale : float, optional
            Uniform scaling factor applied to the output DPI (default ``1``).
        sequence_labels : bool, optional
            If ``True``, draw sequence identifiers to the right of the plot
            (default ``True``).
        legend_entries : list of tuple or None, optional
            Pre-built legend entries as ``(label, color_hex, shape)`` tuples,
            where ``shape`` is ``'rect'``, ``'circle'``, or ``'diamond'``.
            If ``None`` or empty, no legend is drawn (default ``None``).

        Raises
        ------
        ValueError
            If ``sort='tree'`` is requested but no tree was provided at
            initialisation.
        """

        self.plot_type: str = plot_type
        self.output_format: str = output_format

        self.title: str = title
        self._title_font_height: float = self.title_font_size
        self._title_height = 0 if not title else self._title_font_height*2

        self._legend_entries: list = legend_entries or []

        self._ruler_height = 0 if not self.ruler else self._ruler_font_height * 3

        self._plot_floor: float = self.bottom_margin + self._ruler_height

        # Legend sizing — must be calculated before _height so space is reserved
        self._legend_swatch_size: float = self._seq_height
        self._legend_padding: float = self._legend_swatch_size * 0.25
        self._legend_font_size: float = self.seq_name_font_size

        if self._legend_entries:
            self._legend_rows: int = self._compute_legend_rows(self._legend_entries)
            self._legend_height: float = (
                self._legend_rows * (self._legend_swatch_size + self._legend_padding)
                + self._legend_padding
            )
        else:
            self._legend_rows = 0
            self._legend_height = 0

        self._seq_name_width: float = self._max_seq_name_width if sequence_labels else 0
        self._width: float = self.left_margin + self.plot_width + self.plot_label_gap + self._seq_name_width + self.right_margin

        self._height: float = (
            len(self.alignment) * (self._seq_height + self.seq_gap)
            + self.top_margin + self.bottom_margin
            + self._title_height + self._legend_height + self._ruler_height
        )

        self.mark_width: float = mark_width

        self.scale = scale

        self.drawing = Drawing(self._width, self._height)

        if sort == "similar":
            self.sorted_keys = self._sort_similar()
        
        elif sort == "tree":
            if self.tree is None:
                raise ValueError("Cannot sort by tree if no tree is provided")
            
            self.sorted_keys = self._indexes_by_tree_order()

        else: 
            self.sorted_keys = range(len(self.matches_list))

        self._draw_title()
        if self._legend_entries:
            self._draw_legend()
        if self.ruler:
            self._draw_ruler()

        for plot_index, seq_index in enumerate(self.sorted_keys):

            # Add label for sequence
            if sequence_labels:
                id = self.alignment[seq_index].id
                if self.mark_reference:
                    if isinstance(self._mutations.references, int):
                        if seq_index == self._mutations.references:
                            id += " (r)"
                    elif seq_index in self._mutations.references:
                        id += f" (r{self._mutations.references.index(seq_index)+1})"
                    
                x: float = self.left_margin + self.plot_width + self.plot_label_gap #(inch/4)
                y: float = ((self._seq_count-(plot_index + .75)) * (self._seq_height + self.seq_gap))  + self.seq_gap + self._plot_floor
                sequence_str: String = String(x, y, id, fontName="Helvetica", fontSize=self.seq_name_font_size)
                self.drawing.add(sequence_str, id)

            # Add base line for sequence
            color: Color = None

            if plot_type == "match":
                if seq_index in self._mutations.references:
                    color = self._hex_to_color(self._current_scheme[self._mutations.references.index(seq_index)])

                else:
                    sequence: str = ""
                    if isinstance(self.alignment[seq_index], Seq):
                        sequence = str(self.alignment[seq_index])

                    elif isinstance(self.alignment[seq_index], SeqRecord):
                        sequence = str(self.alignment[seq_index].seq)

                    elif isinstance(self.alignment[seq_index], str):
                        sequence = self.alignment[seq_index]

                    sequence = sequence.replace("\n", "")
                    
                    if sequence.replace in self._mutations.references:
                        color = self._hex_to_color(self._current_scheme[self._mutations.references.index(sequence)])

            if not color:
                color = colors.lightgrey

            x1: float = self.left_margin
            x2: float = self.left_margin + self.plot_width
            y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self.seq_gap) + self.seq_gap + self._plot_floor
            sequence_baseline: Line = Line(x1, y, x2, y, strokeColor=color)
            self.drawing.add(sequence_baseline)
    
    def draw_mismatches(self, output_file, *, output_format: str="svg", title: str=None, reference: Union[str, int]=0, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, sort: str="similar", mark_width: float=1, scheme: str="LANL", scale: float=1, sequence_labels: bool=True, legend: bool=True):
        """Generate a mismatch highlighter plot and write it to a file.

        Each sequence track displays coloured marks at positions that differ
        from the reference. Mark colour reflects the residue identity at that
        position according to the chosen colour scheme. Optionally overlays
        symbolic markers for APOBEC mutations (circles), G→A mutations
        (diamonds), stop codons (blue diamonds), and gained/lost
        glycosylation sites (diamonds).

        Parameters
        ----------
        output_file : str or path-like
            Destination file path.
        output_format : str, optional
            Output format string accepted by ``Bio.Graphics._write``, e.g.
            ``'svg'``, ``'pdf'``, ``'png'`` (default ``'svg'``).
        title : str or None, optional
            Title drawn above the plot (default ``None``).
        reference : int or str, optional
            Reference sequence as an alignment index or sequence identifier
            (default ``0``).
        apobec : bool, optional
            Overlay circle markers at APOBEC-signature G→A positions
            (default ``False``).
        g_to_a : bool, optional
            Overlay diamond markers at all G→A positions (default ``False``).
        stop_codons : bool, optional
            Overlay blue diamond markers at stop-codon positions
            (default ``False``).
        glycosylation : bool, optional
            Overlay markers at gained or lost glycosylation sequons
            (default ``False``).
        sort : {'similar', 'tree', 'none'}, optional
            Sequence ordering strategy (default ``'similar'``).
        mark_width : float, optional
            Fractional mark width relative to one alignment column
            (default ``1``).
        scheme : str, optional
            Named colour scheme for residue marks. Built-in options are
            ``'LANL'`` and ``'ML'`` (default ``'LANL'``).
        scale : float, optional
            Output DPI scaling factor (default ``1``).
        sequence_labels : bool, optional
            Draw sequence identifiers to the right of the plot
            (default ``True``).
        legend : bool, optional
            If ``True``, draw a legend above the sequences showing only the
            residues and markers that actually appear in the data
            (default ``True``).

        Returns
        -------
        object
            The return value of ``Bio.Graphics._write`` (format-dependent).
        """

        self._mutations = AlignInfo.Highlighter(self.alignment, seq_type=self.seq_type)

        self.matches_list = self._mutations.list_mismatches(references=reference, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=self.codon_offset)
        self.references = self._mutations.references

        # Scheme must be set before building legend entries and before _setup_drawing
        self.scheme: str = scheme
        self._current_scheme: dict = self.mismatch_plot_colors[self.seq_type][self.scheme] if self.scheme in self.mismatch_plot_colors[self.seq_type] else self.mismatch_plot_colors[self.seq_type]["LANL"]

        self._apobec: bool = apobec
        self._g_to_a: bool = g_to_a
        self._glycosylation: bool = glycosylation
        self._stop_codons: bool = stop_codons

        legend_entries = self._get_legend_entries_mismatch() if legend else []

        self._setup_drawing(output_format=output_format, title=title, sort=sort, mark_width=mark_width, scale=scale, plot_type="mismatch", sequence_labels=sequence_labels, legend_entries=legend_entries)

        for plot_index, seq_index in enumerate(self.sorted_keys):
            matches = self.matches_list[seq_index]
            self._draw_marks_mismatch(plot_index, matches, is_reference=(seq_index == self.references))

        return _write(self.drawing, output_file, self.output_format, dpi=288*self.scale)

    def _draw_marks_mismatch(self, plot_index: int, mismatches: dict[int: list], is_reference: bool=False) -> None:
        """Render mismatch marks and symbolic overlays for a single sequence track.

        Draws filled rectangles for nucleotide/amino acid mismatches and
        overlays symbolic markers (circles or diamonds) for special mutation
        types. Adjacent identical marks are merged into a single wider
        rectangle for visual clarity.

        Parameters
        ----------
        plot_index : int
            Zero-based vertical position of this track in the figure
            (0 = topmost visible track).
        mismatches : dict
            Mismatch annotation dictionary as returned by
            :meth:`Highlighter.get_mismatches_from_str`.
        is_reference : bool, optional
            If ``True``, this track belongs to the reference sequence and
            symbolic glycosylation markers are drawn differently
            (default ``False``).
        """

        already_processed: list = []
        for base, mismatch_item in mismatches.items():
            if base in already_processed:
                continue

            for code in mismatch_item:
                if code in self._current_scheme:
                    width: int = 1
                    color: Color = self._hex_to_color(self._current_scheme[code])

                    check_base: int = base+1
                    while check_base in mismatches and code in mismatches[check_base]:
                        width += 1
                        already_processed.append(check_base)
                        check_base += 1

                    self.drawing.add(self._base_mark(plot_index, base, color, width=width))
        
        # Symboloic markers need to be drawn second so they are on top of the rectangles
        for base, mismatch_item in mismatches.items():
            x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
            y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self.seq_gap) + self.seq_gap + self._plot_floor
                
            if "APOBEC" in mismatch_item:
                self.draw_circle(x, y)
                
            elif "G->A mutation" in mismatch_item:
                self.draw_diamond(x, y)

            elif "Glycosylation" in mismatch_item:
                if is_reference:
                    self.draw_circle(x, y)
                else:
                    if "Glycosylation" not in self.matches_list[self.references].get(base, {}):
                        self.draw_diamond(x, y, filled=True)

            elif "Stop codon" in mismatch_item:
                self.draw_diamond(x, y, color="#0000FF")

        if self.seq_type == "AA" and self._glycosylation:
            for base, mismatch_item in self.matches_list[self.references].items():
                if "Glycosylation" in mismatch_item and "Glycosylation" not in mismatches.get(base, {}):
                    x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
                    y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self.seq_gap) + self.seq_gap + self._plot_floor

                    self.draw_diamond(x, y, color="#0000FF")

    def draw_matches(self, output_file, *, output_format: str="svg", title: str=None, references: list[Union[str, int]]=0, sort: str="similar", mark_width: float=1, scheme: Union[str, dict]="LANL", scale: float=1, sequence_labels: bool=True, legend: bool=True):
        """Generate a match highlighter plot and write it to a file.

        Each sequence track displays coloured marks at positions that are
        identical to one or more reference sequences. Mark colour indicates
        which reference is matched; positions shared with multiple references
        or unique to the query use configurable colours.

        Parameters
        ----------
        output_file : str or path-like
            Destination file path.
        output_format : str, optional
            Output format string (default ``'svg'``).
        title : str or None, optional
            Title drawn above the plot (default ``None``).
        references : int, str, or list thereof, optional
            One or more reference sequences by index or identifier
            (default ``0``).
        sort : {'similar', 'tree', 'none'}, optional
            Sequence ordering strategy (default ``'similar'``).
        mark_width : float, optional
            Fractional mark width relative to one alignment column
            (default ``1``).
        scheme : str or dict, optional
            Colour scheme. Pass a built-in scheme name (``'LANL'`` or
            ``'ML'``) or a dictionary with keys ``'references'``
            (list of hex colour strings), ``'unique'`` (hex string or
            ``None``), and ``'multiple'`` (hex string or ``None``).
            Default is ``'LANL'``.
        scale : float, optional
            Output DPI scaling factor (default ``1``).
        sequence_labels : bool, optional
            Draw sequence identifiers to the right of the plot
            (default ``True``).
        legend : bool, optional
            If ``True``, draw a legend above the sequences showing only the
            residues and markers that actually appear in the data
            (default ``True``).

        Returns
        -------
        object
            The return value of ``Bio.Graphics._write``.

        Raises
        ------
        ValueError
            If a string scheme name is not recognised, or if a dict scheme
            is missing required keys or empty reference colours.
        TypeError
            If ``scheme`` is neither a str nor a dict.
        """

        self._mutations = AlignInfo.Highlighter(self.alignment, seq_type=self.seq_type)

        self.matches_list = self._mutations.list_matches(references=references)
        self.references = self._mutations.references

        # Scheme must be set before building legend entries and before _setup_drawing
        if isinstance(scheme, str):
            self.scheme: str = scheme

            if self.scheme not in self.match_plot_colors:
                raise ValueError(f"Scheme {self.scheme} is not a valid scheme")

            self._current_scheme: list = self.match_plot_colors[self.scheme]["references"]
            self._current_unique_color: str = self.match_plot_colors[self.scheme]["unique"]
            self._current_multiple_color: str = self.match_plot_colors[self.scheme]["multiple"]

        elif isinstance(scheme, dict):
            if "references" not in scheme or "unique" not in scheme or "multiple" not in scheme:
                raise ValueError("Scheme dictionary must contain 'references', 'unique', and 'multiple' keys")

            if not scheme["references"]:
                raise ValueError("Scheme dictionary must contain at least one 'reference' color")

            if "multiple" not in scheme or (scheme["multiple"] is not None and not scheme["multiple"]):
                raise ValueError("Scheme dictionary must contain a 'multiple' color")

            if "unique" not in scheme or (scheme["unique"] is not None and not scheme["unique"]):
                raise ValueError("Scheme dictionary must contain a 'unique' color")

            self._current_scheme: list = scheme["references"]
            self._current_unique_color: str = scheme["unique"]
            self._current_multiple_color: str = scheme["multiple"]

        else:
            raise TypeError("scheme must be a string or a dictionary")

        legend_entries = self._get_legend_entries_match() if legend else []

        self._setup_drawing(output_format=output_format, title=title, sort=sort, mark_width=mark_width, scale=scale, plot_type="match", sequence_labels=sequence_labels, legend_entries=legend_entries)

        for plot_index, seq_index in enumerate(self.sorted_keys):
            matches = self.matches_list[seq_index]
            self._draw_marks_match(plot_index, matches, is_reference=(seq_index in self.references))

        return _write(self.drawing, output_file, self.output_format, dpi=288*self.scale)

    def _get_legend_entries_mismatch(self) -> list[tuple]:
        """Build ordered legend entries from the mismatch data actually present.

        Scans ``self.matches_list`` to find which residues and special
        annotation types appear in the data. For nucleotide alignments
        residues are listed in canonical order (A, C, G, T/U, Gap). For
        amino acid alignments, residues sharing the same colour in the active
        scheme are grouped into a single entry with a combined label
        (e.g. ``'I/L/V'``), greatly reducing legend clutter.

        Symbolic marker entries (APOBEC, G→A, stop codon, glycosylation) are
        appended after residue entries, but only when the corresponding flag
        was enabled and the marker actually appears in the data.

        Returns
        -------
        list of tuple
            Ordered ``(label, color_hex, shape)`` tuples where ``shape`` is
            one of ``'rect'``, ``'circle'``, or ``'diamond'``.
        """

        # Collect residues and special codes that actually appear
        seen_residues: set = set()
        seen_specials: set = set()

        for seq_mismatches in self.matches_list:
            for annotations in seq_mismatches.values():
                seen_residues.add(annotations[0])
                for code in annotations[1:]:
                    seen_specials.add(code)

        entries: list = []

        if self.seq_type == "NT":
            nt_order = ["A", "C", "G", "T", "U", "Gap"]
            for residue in nt_order:
                if residue in seen_residues and residue in self._current_scheme:
                    entries.append((residue, self._current_scheme[residue], "rect"))

        else:  # AA — group residues by colour
            # Invert the colour map: hex -> [residues with that colour]
            color_to_residues: dict = {}
            for residue, hex_color in self._current_scheme.items():
                if residue in seen_residues:
                    color_to_residues.setdefault(hex_color, []).append(residue)

            # Preserve scheme insertion order for group ordering
            seen_colors: list = []
            for residue in self._current_scheme:
                hex_color = self._current_scheme[residue]
                if hex_color in color_to_residues and hex_color not in seen_colors:
                    seen_colors.append(hex_color)

            for hex_color in seen_colors:
                residues = color_to_residues[hex_color]
                label = "/".join(sorted(residues))
                entries.append((label, hex_color, "rect"))

            if "Gap" in seen_residues and "Gap" in self._current_scheme:
                entries.append(("Gap", self._current_scheme["Gap"], "rect"))

        # Symbolic markers — only if enabled and present in data
        if self._g_to_a and "G->A mutation" in seen_specials:
            entries.append(("G→A", "#FF00FF", "diamond"))

        if self._apobec and "APOBEC" in seen_specials:
            entries.append(("APOBEC", "#FF00FF", "circle"))

        if self._stop_codons and "Stop codon" in seen_specials:
            entries.append(("Stop codon", "#0000FF", "diamond"))

        if self._glycosylation and "Glycosylation" in seen_specials:
            entries.append(("Glycosylation (gained)", "#FF00FF", "diamond_filled"))
            entries.append(("Glycosylation (ref)", "#FF00FF", "circle"))
            entries.append(("Glycosylation (lost)", "#0000FF", "diamond"))

        return entries

    def _get_legend_entries_match(self) -> list[tuple]:
        """Build ordered legend entries from the match data actually present.

        Scans ``self.matches_list`` to find which reference indices, and
        whether ``'Unique'`` or ``'Multiple'`` categories appear. Entries are
        ordered R1, R2, … followed by Multiple then Unique, omitting any that
        do not appear in the data.

        Returns
        -------
        list of tuple
            Ordered ``(label, color_hex, shape)`` tuples.
        """

        seen_refs: set = set()
        has_unique: bool = False
        has_multiple: bool = False

        for seq_matches in self.matches_list:
            for annotations in seq_matches.values():
                if "Unique" in annotations:
                    has_unique = True
                elif len(annotations) > 1:
                    has_multiple = True
                else:
                    seen_refs.update(i for i in annotations if isinstance(i, int))

        entries: list = []

        for ref_index in sorted(seen_refs):
            if ref_index < len(self._current_scheme) and self._current_scheme[ref_index] is not None:
                label = f"R{ref_index + 1}"
                entries.append((label, self._current_scheme[ref_index], "rect"))

        if has_multiple and self._current_multiple_color is not None:
            entries.append(("Multiple", self._current_multiple_color, "rect"))

        if has_unique and self._current_unique_color is not None:
            entries.append(("Unique", self._current_unique_color, "rect"))

        return entries

    def _compute_legend_rows(self, entries: list) -> int:
        """Return the number of rows needed to lay out legend entries across the plot width.

        Each entry occupies a swatch plus its label text. Entries are packed
        left-to-right and wrap onto a new row when they would exceed
        ``self.plot_width``.

        Parameters
        ----------
        entries : list of tuple
            Legend entries as ``(label, color_hex, shape)`` tuples.

        Returns
        -------
        int
            Number of rows required (minimum 1 if entries is non-empty).
        """

        x: float = 0
        rows: int = 1
        swatch: float = self._legend_swatch_size
        gap: float = self._legend_padding

        for label, _, _ in entries:
            label_width: float = stringWidth(label, self.seq_name_font, self._legend_font_size)
            entry_width: float = swatch + gap + label_width + gap * 2

            if x + entry_width > self.plot_width and x > 0:
                rows += 1
                x = 0

            x += entry_width

        return rows

    def _draw_legend(self) -> None:
        """Render the legend between the title and the sequence tracks.

        Lays out coloured swatches (rectangles, circles, or diamonds) with
        text labels in rows across the plot width. The legend is positioned
        immediately below the title (or below the top margin if there is no
        title), and the sequence tracks begin below it.

        The vertical position of the legend is derived from ``_height``,
        ``top_margin``, and ``_title_height``, matching the coordinate system
        used by ``_draw_title`` and the sequence baselines.
        """

        swatch: float = self._legend_swatch_size
        pad: float = self._legend_padding
        font_size: float = self._legend_font_size

        # Top of the legend block in drawing coordinates
        legend_top: float = self._height - self.top_margin - self._title_height - pad

        x: float = self.left_margin
        row: int = 0

        for label, hex_color, shape in self._legend_entries:
            label_width: float = stringWidth(label, self.seq_name_font, font_size)
            entry_width: float = swatch + pad + label_width + pad * 2

            if x + entry_width > self.left_margin + self.plot_width and x > self.left_margin:
                row += 1
                x = self.left_margin

            # Vertical centre of this row's swatch
            y_centre: float = legend_top - (row * (swatch + pad)) - swatch / 2

            color: Color = self._hex_to_color(hex_color)

            if shape == "rect":
                self.drawing.add(Rect(
                    x, y_centre - swatch / 2,
                    swatch, swatch,
                    fillColor=color, strokeColor=color, strokeWidth=0.1,
                ))
            elif shape == "circle":
                self.drawing.add(Circle(
                    x + swatch / 2, y_centre, swatch / 2,
                    fillColor=color, strokeColor=color, strokeWidth=0.1,
                ))
            elif shape == "diamond":
                half: float = swatch / 2
                self.drawing.add(Polygon(
                    [x + half, y_centre - half,
                     x,        y_centre,
                     x + half, y_centre + half,
                     x + swatch, y_centre],
                    strokeColor=color, strokeWidth=2, fillColor=None,
                ))
            elif shape == "diamond_filled":
                half: float = swatch / 2
                self.drawing.add(Polygon(
                    [x + half, y_centre - half,
                     x,        y_centre,
                     x + half, y_centre + half,
                     x + swatch, y_centre],
                    strokeColor=color, strokeWidth=2, fillColor=color,
                ))

            # Label to the right of the swatch
            self.drawing.add(String(
                x + swatch + pad, y_centre - font_size / 2,
                label,
                fontName=self.seq_name_font, fontSize=font_size,
            ))

            x += entry_width

    def _draw_marks_match(self, plot_index: int, matches: dict[int, list], is_reference: bool) -> None:
        """Render match marks for a single sequence track.

        Draws coloured rectangles at positions where the query matches a
        reference. Adjacent identical marks are merged into a single wider
        rectangle. Three mutually exclusive mark categories are handled in
        priority order: unique positions, multi-reference matches, and
        single-reference matches.

        Parameters
        ----------
        plot_index : int
            Zero-based vertical position of this track in the figure.
        matches : dict
            Match annotation dictionary as returned by
            :meth:`Highlighter.get_matches_from_str`.
        is_reference : bool
            Unused currently; reserved for future reference-track styling.
        """

        already_processed: dict[list] = {"Unique": [], "Multiple": [], "Single": []}
        for base, match_item in matches.items():
            check_base: int = base+1
            width: int = 1

            if "Unique" in match_item:
                if base in already_processed["Unique"]:
                    continue

                if self._current_unique_color is not None:
                    color: Color = self._hex_to_color(self._current_unique_color)

                    while check_base in matches and "Unique" in matches[check_base]:
                        width += 1
                        already_processed["Unique"].append(check_base)
                        check_base += 1

                    self.drawing.add(self._base_mark(plot_index, base, color, width=width))

            elif len(match_item) > 1:

                if base in already_processed["Multiple"]:
                    continue

                if self._current_multiple_color is not None:
                    color: Color = self._hex_to_color(self._current_multiple_color)

                    while check_base in matches and len(matches[check_base]) > 1:
                        width += 1
                        already_processed["Multiple"].append(check_base)
                        check_base += 1

                    self.drawing.add(self._base_mark(plot_index, base, color, width=width))

            else:
                if base in already_processed["Single"]:
                    continue

                for code in match_item:
                    if self._current_scheme[code] is not None:
                        color: Color = self._hex_to_color(self._current_scheme[code])

                        while check_base in matches and code in matches[check_base]:
                            width += 1
                            already_processed["Single"].append(check_base)
                            check_base += 1

                        self.drawing.add(self._base_mark(plot_index, base, color, width=width))

    def _base_mark(self, plot_index: int, base: int, color: Color, width: int=1) -> Rect:
        """Construct a filled rectangle mark for one or more consecutive bases.

        Parameters
        ----------
        plot_index : int
            Zero-based vertical position of this track in the figure.
        base : int
            Zero-based alignment position of the leftmost column to mark.
        color : Color
            ReportLab ``Color`` instance for both the fill and stroke.
        width : int, optional
            Number of consecutive alignment columns to span (default ``1``).

        Returns
        -------
        Rect
            A ReportLab ``Rect`` object ready to be added to the drawing.
        """

        x1: float = self.left_margin + self._base_left(base)
        x2: float = self.left_margin + self._base_left(base + (width - 1) + self.mark_width )

        y1: float = ((self._seq_count-plot_index) * (self._seq_height + self.seq_gap)) + (self.seq_gap/2) + self._plot_floor
        y2: float = ((self._seq_count-(plot_index+1)) * (self._seq_height + self.seq_gap)) + self.seq_gap + self._plot_floor

        return Rect(x1, y1, x2-x1, y2-y1, fillColor=color, strokeColor=color, strokeWidth=0.1)

    def _draw_title(self) -> None:
        """Draw the plot title centred above the sequence tracks.

        Has no effect if ``self.title`` is ``None`` or an empty string.
        The title is positioned relative to the top margin using the
        configured title font and font size.
        """
            
        if self.title:
            x: float = self.left_margin + (self.plot_width/2)
            y: float = self._height - self.top_margin - (self._title_font_height/2)

            self.drawing.add(String(x, y, self.title, textAnchor="middle", fontName=self.title_font, fontSize=self.title_font_size))

    def _draw_ruler(self) -> None:
        """Draw the position ruler along the bottom of the plot.

        For alignments of 20 columns or fewer every position is labelled.
        For longer alignments, major tick positions are distributed evenly
        across the sequence, rounded to significant digits, with minor ticks
        interpolated between them. A solid horizontal line is drawn across
        the full plot width at the top of the tick marks.
        """

        if self._seq_length <= 20:
            self._ruler_marks(range(self._seq_length))
        else:
            bases = [0]
            last_base: int = self.significant_digits(self._seq_length)-1
            
            dist = last_base / self.ruler_major_ticks

            for base in range(1, self.ruler_major_ticks+1):
                bases.append(int(base * dist))

            bases.append(last_base)

            self._ruler_marks(bases)

        # Draw vertical line
        x1: float = self.left_margin
        x2: float = self.left_margin + self.plot_width
        y: float = self.bottom_margin + (self._ruler_font_height * 3)

        ruler_line: Line = Line(x1, y, x2, y, strokeColor=colors.black, strokeWidth=2)
        self.drawing.add(ruler_line)

    def _ruler_marks(self, marked_bases: list) -> None:
        """Draw major tick marks and labels, and interpolate minor ticks between them.

        Parameters
        ----------
        marked_bases : list of int
            Zero-based alignment positions at which to place major ticks
            and numeric labels.
        """

        mark_locations: list = []

        for index, base in enumerate(marked_bases):
                mark_locations.append(self._ruler_label(base))
                self._ruler_heavy_tick(base)

                if index:
                    spacing = (mark_locations[index][0] - mark_locations[index-1][0]) / (self.ruler_minor_ticks+1)
                    left = mark_locations[index-1][0]

                    for tick in range(1, self.ruler_minor_ticks+1):
                        self._ruler_light_tick(left + (tick * spacing))
                    

    def _ruler_label(self, base: int) -> tuple[float, float]:
        """Draw a 1-based numeric label on the ruler at the given alignment position.

        Parameters
        ----------
        base : int
            Zero-based alignment position. The displayed label is ``base + 1``.

        Returns
        -------
        tuple of float
            ``(x, y)`` coordinates of the label anchor in ReportLab units.
        """

        x: float = self.left_margin + self._base_center(base)
        y: float = self.bottom_margin+self._ruler_font_height

        self.drawing.add(String(x, y, str(base + 1), textAnchor="middle", fontName=self.ruler_font, fontSize=self.ruler_font_size))

        return (x, y)

    def _ruler_heavy_tick(self, base: int) -> tuple[float, float, float]:
        """Draw a major (heavy) tick mark on the ruler at the given alignment position.

        Parameters
        ----------
        base : int
            Zero-based alignment position.

        Returns
        -------
        tuple of float
            ``(x, top, bottom)`` coordinates of the tick line in ReportLab units.
        """

        x: float = self.left_margin + self._base_center(base)
        top: float = self.bottom_margin + (self._ruler_font_height * 3)
        bottom: float = self.bottom_margin + (self._ruler_font_height * 2)

        self.drawing.add(Line(x, top, x, bottom, strokeColor=colors.black, strokeWidth=1))

        return (x, top, bottom)

    def _ruler_light_tick(self, x: float) -> tuple[float, float, float]:
        """Draw a minor (light) tick mark on the ruler at the given x coordinate.

        Unlike :meth:`_ruler_heavy_tick`, this method accepts an absolute
        x coordinate (including the left margin) rather than an alignment
        position, because minor ticks are interpolated between major tick
        positions.

        Parameters
        ----------
        x : float
            Absolute x coordinate in ReportLab units (left margin already
            included).

        Returns
        -------
        tuple of float
            ``(x, top, bottom)`` coordinates of the tick line in ReportLab units.
        """

        top: float = self.bottom_margin + (self._ruler_font_height * 3)
        bottom: float = self.bottom_margin + (self._ruler_font_height * 2.5)

        self.drawing.add(Line(x, top, x, bottom, strokeColor=colors.black, strokeWidth=1))

        return (x, top, bottom)

    def _base_left(self, base: int) -> float:
        """Return the left edge x coordinate of an alignment column, relative to the plot origin.

        The value does *not* include ``left_margin``; callers must add it
        explicitly when placing elements on the drawing.

        Parameters
        ----------
        base : int
            Zero-based alignment column index.

        Returns
        -------
        float
            X coordinate in ReportLab units, measured from the left edge of
            the plot area (i.e. excluding the left margin).
        """

        # if base == self._seq_length:
        #     return self.plot_width
        
        return (base / self._seq_length) * self.plot_width
    
    def _base_center(self, base: int) -> float:
        """Return the centre x coordinate of an alignment column, relative to the plot origin.

        Parameters
        ----------
        base : int
            Zero-based alignment column index.

        Returns
        -------
        float
            X coordinate in ReportLab units, measured from the left edge of
            the plot area (i.e. excluding the left margin).
        """

        left_x: float = self._base_left(base)
        right_x: float = self._base_left(base + 1)

        return left_x + ((right_x-left_x)/2)
    
    @property
    def _max_seq_name_width(self) -> float:
        """Return the rendered pixel width of the longest sequence label.

        Accounts for the reference indicator suffix (``' (r)'`` for mismatch
        plots; ``' (r10)'`` for match plots, sizing for up to 10 references)
        so that the label column is wide enough for all sequences including
        the reference.

        Returns
        -------
        float
            Maximum label width in ReportLab units as reported by
            ``reportlab.pdfbase.pdfmetrics.stringWidth``.
        """

        max_width: float = 0
        if self.plot_type == "match":
            reference_tag: str = " (r10)"
        elif self.plot_type == "mismatch":
            reference_tag: str = " (r)"

        for sequence in self.alignment:
            width: float = stringWidth(f"{sequence.id} {reference_tag}", self.seq_name_font, self.seq_name_font_size)
            if width > max_width:
                max_width = width
        
        return max_width
    
    def _hex_to_color(self, hex: str) -> Color:
        """Convert a CSS hex colour string to a ReportLab ``Color`` object.

        Parameters
        ----------
        hex : str
            Six-digit hexadecimal colour string, with or without a leading
            ``'#'`` (e.g. ``'#FF0000'`` or ``'FF0000'``).

        Returns
        -------
        Color
            ReportLab ``Color`` with ``r``, ``g``, ``b`` components in the
            range [0, 1].
        """

        hex = hex.lstrip("#")
        color_list = [int(hex[i:i+2], 16)/256 for i in (0, 2, 4)]
        return Color(color_list[0], color_list[1], color_list[2])
    
    def _sort_similar(self) -> list[int]:
        """Return alignment indices sorted by ascending number of annotated positions.

        Sequences with fewer annotated positions (i.e. more similar to the
        reference) are placed first. This mirrors the display convention of
        the LANL highlighter tool.

        Returns
        -------
        list of int
            Alignment indices in ascending order of ``len(matches_list[i])``.
        """

        return sorted(range(len(self.matches_list)), key=lambda x: len(self.matches_list[x]))
    
    def draw_diamond(self, x: float, y: float, color: str="#FF00FF", filled: bool=False) -> None:
        """Draw a diamond-shaped marker on the plot at the given coordinates.

        Used to overlay symbolic annotations such as G→A mutations, gained
        glycosylation sites, and stop codons on top of sequence tracks.

        Parameters
        ----------
        x : float
            Absolute x coordinate of the diamond centre in ReportLab units.
        y : float
            Absolute y coordinate of the diamond centre in ReportLab units.
        color : str, optional
            Hex colour string for the diamond stroke (and fill, if
            ``filled=True``) (default ``'#FF00FF'``).
        filled : bool, optional
            If ``True``, the diamond interior is filled with ``color``;
            otherwise it is transparent (default ``False``).
        """
        
        fill_color = self._hex_to_color(color) if filled else None

        diamond = Polygon([x, y-((self._seq_height/3)/2), x-((self._seq_height/3)/2), y, x, y+((self._seq_height/3)/2), x+((self._seq_height/3)/2), y], strokeColor=self._hex_to_color(color), strokeWidth=2, fillColor=fill_color)
        
        self.drawing.add(diamond)

    def draw_circle(self, x: float, y: float, color: str="#FF00FF", filled: bool=True) -> None:
        """Draw a circle marker on the plot at the given coordinates.

        Used to overlay symbolic annotations such as APOBEC mutations and
        reference glycosylation sites on top of sequence tracks.

        Parameters
        ----------
        x : float
            Absolute x coordinate of the circle centre in ReportLab units.
        y : float
            Absolute y coordinate of the circle centre in ReportLab units.
        color : str, optional
            Hex colour string for the circle fill (default ``'#FF00FF'``).
        filled : bool, optional
            Currently unused; the circle is always drawn filled
            (default ``True``).
        """
        
        circle = Circle(x, y, (self._seq_height/3)/2, fillColor=self._hex_to_color(color), strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=0.1)
        self.drawing.add(circle)
    
    def _get_index_by_id(self, id: str) -> int:
        """Return the alignment index of a sequence by its identifier.

        Parameters
        ----------
        id : str
            The sequence identifier to search for (matches ``SeqRecord.id``).

        Returns
        -------
        int
            Zero-based index of the matching sequence in the alignment.

        Raises
        ------
        IndexError
            If no sequence with the given identifier is found.
        """

        for index, sequence in enumerate(self.alignment):
            if sequence.id == id:
                return index
        
        raise IndexError(f"Sequence with id '{id}' not found")
    
    def _indexes_by_tree_order(self) -> list[int]:
        """Return alignment indices in the order of terminal nodes in the phylogenetic tree.

        Used when ``sort='tree'`` is requested in a draw method. Terminal
        names in the tree must match sequence identifiers in the alignment.

        Returns
        -------
        list of int
            Alignment indices ordered according to ``self.tree.get_terminals()``.

        Raises
        ------
        IndexError
            If a terminal name in the tree cannot be found in the alignment.
        """

        indexes: list[int] = []

        for terminal in self.tree.get_terminals():
            indexes.append(self._get_index_by_id(terminal.name))

        return indexes
    
    @staticmethod
    def significant_digits(value: float, digits: int=2) -> float:
        """Round a positive integer down to a given number of significant digits.

        Used to produce round ruler tick positions (e.g. 987 → 980 with two
        significant digits). Only positive values are rounded; zero and
        negative values are returned unchanged.

        Parameters
        ----------
        value : float
            The value to round.
        digits : int, optional
            Number of significant digits to retain (default ``2``).

        Returns
        -------
        float
            The rounded value, or ``value`` unchanged if it has fewer digits
            than ``digits`` or is non-positive.

        Examples
        --------
        >>> HighlighterPlot.significant_digits(987)
        980
        >>> HighlighterPlot.significant_digits(1234, digits=3)
        1230
        """

        if value == 0:
            return 0

        if value > 0:
            value_str = str(value)
            value_len = len(value_str)

            if value_len > digits:
                return(int(value_str[:digits] + "0" * (value_len-digits)))
            else:
                return value

        return value
    
    @staticmethod
    def guess_alignment_type(alignment) -> str:
        """Infer whether an alignment contains nucleotide or amino acid sequences.

        Scans every sequence for characters that fall outside the IUPAC
        nucleotide alphabet. If any such character is found the alignment is
        classified as amino acid; otherwise it is classified as nucleotide.

        Parameters
        ----------
        alignment : list-like of str, Seq, or SeqRecord
            The sequences to inspect.

        Returns
        -------
        str
            ``'NT'`` if all characters are valid IUPAC nucleotide codes,
            ``'AA'`` otherwise.

        Notes
        -----
        The IUPAC nucleotide codes recognised are: ``A C G T U i R Y K M S W
        B D H V N -``. The gap character ``'-'`` is treated as a nucleotide
        code so that pre-aligned sequences are handled correctly.
        """

        nt_codes: str = "ACGTUiRYKMSWBDHVN-" # IUPAC nucleotide codes

        for sequence in alignment:
            if isinstance(sequence, SeqRecord):
                sequence_str = str(sequence.seq)
            elif isinstance(sequence, Seq):
                sequence_str = str(sequence)
            elif isinstance(sequence, str):
                sequence_str = sequence
                
            for symbol in sequence_str.strip():
                if symbol not in nt_codes:
                    return "AA"
                    
        return "NT"

Graphics.HighlighterPlot = HighlighterPlot


from Bio import SeqUtils

@cache
def codon_position(sequence: Union[str, Seq, SeqRecord], base: int, *, codon_offset: int=0) -> int:
    """Return the codon position of a residue in an alignment, accounting for gaps.

    Determines whether the residue at ``base`` is the 1st, 2nd, or 3rd
    position within its codon. Gap characters (``'-'``) in the sequence
    before ``base`` are excluded from the position count so that the reading
    frame is preserved across indels.

    This function is registered into the Biopython namespace as
    ``Bio.SeqUtils.codon_position`` and is also used internally by
    :meth:`Highlighter.get_mismatches_from_str` for stop codon detection.

    Results are cached via :func:`functools.cache` for performance.

    Parameters
    ----------
    sequence : str, Seq, or SeqRecord
        The aligned nucleotide sequence containing the residue of interest.
    base : int
        Zero-based index of the residue within ``sequence``.
    codon_offset : int, optional
        Reading frame offset to apply (default ``0``). A value of ``1``
        indicates the first base of the sequence is the second position of
        a codon, and so on. The offset is applied after gap adjustment.

    Returns
    -------
    int
        Codon position: ``0`` for the first base of a codon, ``1`` for the
        second, or ``2`` for the third.

    Raises
    ------
    TypeError
        If ``sequence`` is not a str, Seq, or SeqRecord, or if ``base`` is
        not an int.
    ValueError
        If ``base`` exceeds the length of ``sequence``, or if the residue at
        ``base`` is a gap character.

    Examples
    --------
    >>> from Bio import SeqUtils
    >>> SeqUtils.codon_position("ATG---TTA", 6)
    0
    >>> SeqUtils.codon_position("ATG---TTA", 7)
    1
    """

    if isinstance(sequence, Seq):
            sequence = str(sequence)
    elif isinstance(sequence, SeqRecord):
            sequence = str(sequence.seq)
    elif not isinstance(sequence, str):
            raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")
    
    if not isinstance(base, int):
        raise TypeError(f"Expected position to be an int, got {type(base)}")

    if base > len(sequence):
        raise ValueError(f"Position {base} is greater than the length of the sequence ({len(sequence)})")
    
    if sequence[base] == "-":
        raise ValueError(f"Position {base} is a gap")
    
    adjusted_base: int =  base-sequence[:base+1].count("-")
    return ((adjusted_base + codon_offset) % 3)

SeqUtils.codon_position = codon_position