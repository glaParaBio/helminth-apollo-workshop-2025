# Exploring annotation 

In an Apollo session open the *Trichuris trichiura* assembly. Click on "open
track selector" and select the "Annotations" and "T. trichiura IsoSeq" tracks.

Using the coordinate search box navigate to position
`TTRE_chr2:10,913,100-10,916,800`. Explore this regions by moving left and
right, zoom in and out. Then return to this position.

You may notice that the 3'UTR of gene TTRE4940 extends too far relative to the
IsoSeq data so let's fix it by dragging the leftmost edge of the exon to match
the IsoSeq reads. We should also resize the underlying transcript and gene. For
this, right-click on an intron of this gene and select "Edit feature details".
Edit the "End" of this mRNA to match the end of the exon you just edited. You
can expand the "Related features" box to get the exon coordinates. The select
the parent feature (i.e. the gene) and again edit its end coordinate to match
the transcript. We may want to add an attribute to this transcript to describe
what we did.

Next search for "TTRE688" using the coordinate box. This time the transcript
needs a UTR to be extended. We can use the same procedure as above, but this
time we need first to extend the gene, then transcript.

Now search for TTRE6744. Since there is weak support from IsoSeq we may
consider deleting this gene altogether. We can use the "Delete feature"
operation accessible from right-clicking on the transcript. Alternatively, we
can display the table editor and do it from there. To open the table editor,
click on the three vertical dots next to the track name on the left of the
screen. Then from sub-menu Appearance select "Show both graphical and Table
display". Right-click on the row of this gene and delete it.

We can also create a new feature from scratch.  Click and drag an area in the
ruler and then select "Add new feature". Look at the possible ontology types
suggested, and then choose "match". On the new feature, click "Add child
feature" and choose the coordinates and "Match part" as the type. Now
right-click on the "match" feature and select "Delete feature".

NB: If you want to explore these data by yourself note that to keep the data
files small we only included genes and synteny links in these intervals:

*T. muris*: TMUE_LG2:17202488-20000820

*T. trichiura*: TTRE_chr2:10695245-13020708

# Using synteny

We will now make use of synteny to help deciding what genes should be edited
and how. 

Start a new session by selecting "New session" from the "File" menu.

Select "Linear synteny view" from the dropdown menu. In row 1 select *Trichuris muris*
and in row 2 *Trichuris trichiura*. Now the radio button "Existing track"
should be checked and the available track should be
"trichuris_muris_vs_trichuris_trichiura TBLASTX" we loaded earlier. Click "Launch".
The first view is fully zoomed out and saturated with synteny links. 

In the coordinate box on the right, the one corresponding to *T.  trichiura*,
search for TTRE2050 and zoom out a bit. Now let's open the IsoSeq track for
*trichiura*. In synteny view we can access the menu for each track by clicking
on the three vertical dots on the left, next to the search box, then select
"Row view menus", "View menu 2", "Open track selector". Also zoom in on the
*muris* track to inspect the genes in the region syntenic to TTRE2050. We may
decide that Transcript TTRE2050_t2 has weak support from IsoSeq or synteny and
we may want to edit or delete it.

# Regions of interest

These are some example genes in *T. trichiura* that probably need some curation:

```
TTRE4940 Shorten UTR *
TTRE688 Extend UTR *
TTRE6744 Delete *
TTRE6731 Delete gene *
TTRE6830 Delete gene *
TTRE2050 Delete tx TMUE_LG2:18,234,381..18,272,122[rev] TTRE_chr2:11,715,666..11,730,260
TTRE6749 Add UTR exon also CDS TMUE_LG2:17,211,180..17,216,769[rev] TTRE_chr2:10,946,675..10,951,590
TTRE6722 TMUE_LG2:18,046,592..18,050,363 TTRE_chr2:10,721,018..10,724,789
         TMUE_LG2:17,435,215..17,440,152[rev]
TTRE6839 Shorten CDS TMUE_LG2:19,449,863..19,458,308[rev] TTRE_chr2:12,163,993..12,172,437
TTRE6847 Add UTR (possibly also CDSs) TMUE_LG2:19,528,245..19,531,940 TTRE_chr2:12,227,640..12,238,090
```
