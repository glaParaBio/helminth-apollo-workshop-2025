# Exploring annotation 

Navigate your browser to [http://localhost/](http://localhost/) and select
"Linear Genome View". Then select the *Trichuris trichiura* assembly from the
dropdown menu and open it.

Let's visualize the annotation track to be curated and the IsoSeq track. Click
on "Open track selector" and check the boxes for "Annotations" and "T.
trichiura IsoSeq" tracks.

Using the coordinate search box navigate to position
`TTRE_chr2:10,913,100-10,916,800`. Explore this regions by moving left and
right, zooming in and out, dragging the track left or right. Then return to
this position.

## Making a UTR shorter

You may notice that the 3'UTR of gene TTRE4940
(TTRE_chr2:10,913,100-10,916,800) extends too far relative to the IsoSeq data.
Let's fix it by dragging the leftmost edge of the exon to match the IsoSeq
reads (you may want to zoom into the 3'UTR for a better view). 

We should also resize the underlying transcript and gene to match the new end
of the UTR. For this, right-click on an intron of this transcript and select
"Edit feature details".  Edit the "End" of this mRNA to match the end of the
exon you just edited (you can expand the "Related features" box to see the exon
coordinates).  From the expand "Related features" box select the parent feature
(i.e. the gene) and again edit its end coordinate to match the transcript. 

Once we are satisfied with the edits, we may want to add a new attribute to
this transcript to describe what we did. For this use the "Add new" option
under the Attributes section.

## Making a UTR longer

Next search for "TTRE688" using the coordinate box. Zoom into the 5'UTR (note
this gene is on the reverse strand) and you will notice that the UTR to be
extended. Pop up the "Edit feature details" as above and use the same procedure
to correct this gene. However, this time we need to first extend the gene, then
the transcript since we cannot have a child feature extending beyond its parent
(*e.g.* an exon cannot extend beyond its parent transcript. The Apollo
developers are reconsidering this behaviour).

## Deleting a gene

Now search for TTRE6744. Since there is weak support from IsoSeq we may delete
this gene altogether. We can use the "Delete feature" operation accessible from
right-clicking on the transcript. Alternatively, we can display the table
editor and do it from there. To open the table editor, click on the three
vertical dots next to the track name on the left of the screen. From sub-menu
"Appearance" select "Show both graphical and Table display". You can make edits
in this table to correct coordinates, feature types, or in this case delete
features. Right-click on the row of this gene and delete it.

## Creating a feature from scratch

We can also create a new feature from scratch.  Click and drag an area in the
ruler and then select "Add new feature". Look at the possible ontology types
suggested, and then choose "match". On the new feature, click "Add child
feature" and choose the coordinates and "Match part" as the type. Now
right-click on the "match" feature and select "Delete feature".

# Using synteny

We will now make use of synteny to help deciding what genes should be edited
and how.

Start a new session by selecting "New session" from the "File" menu.

Select "Linear synteny view" from the dropdown menu. In row 1 select *Trichuris muris*
and in row 2 *Trichuris trichiura*. Now the radio button "Existing track"
should be checked and the available track should be
"trichuris_muris_vs_trichuris_trichiura TBLASTX" that we loaded earlier. Click "Launch".
The first view is fully zoomed out and saturated with synteny links. Click on
"Open track selector" for the top track and check Annotation, do the same for
the bottom track and also check the IsoSeq track.

First let's explore some functionalities of the synteny view. In the *T muris*
search box (the one on the left) go to region TMUE_LG2:18,440,431..18,443,422
and you should see a few synteny links. Right-click on one of them and click on
"Center on feature" to center both track on the region covered by this link. 

Since this view may be too narrow you can zoom out each track inidividually
(see lenses on the right) or alternatively you can "Link views". To link views
click on the the three vertical dots on the left of the search boxes and check
the option "Link views". Now the zooming will apply to both tracks. 

Note that you clicking and dragging on the synteny track (the one in the
middle) will move both genome tracks together even if views are unlinked.

After zooming out a few times you will notice that region (approximately)
TMUE_LG2:18,439,000..18,444,000 (*muris*) is inverted relative to the
corresponding TTRE_chr2:11,457,000..11,462,000 (*trichiura*). We can flip one
of the views, say the first, by clicking again the three vertical dots, then
"Row views menu", "View menu 1", "Horizontally flip". Now the synteny picture
should be clearer. Remember to unlink views if you want to move one track but
not the other.

Once you are here, you may consider a small adjustment to the UTRs of TTRE6784
to match the IsoSeq data.

Now search for the *trichiura* gene TTRE6839 and adjust the overall view as you
see fit. The long CDS on the 3'end of the transcript doesn't seem very
convincing based on the evidence and we should delete it. These are the steps I
would suggest, but you can do the same in various ways:

* Zoom in on the preceeding CDS, approximately TTRE_chr2:12,165,988..12,166,272
* Extend this CDS to a suitable stop codon using the IsoSeq as guide (zoom in more as needed)
* (You should see a warning indicating that presence of an internal stop codon, this is expected)
* Right click on the CDS and bring up the "Edit feature details" dialog. Edit
  the end of CDS to match the position of the stop codon you selected above. (The warning should disappear)
* Delete the 3' exon that now looks like a UTR
* As before, resize transcript and gene accordingly
* Add an attribute to describe the change

<!---
# Regions of interest

These are some example genes in *T. trichiura* that probably need some curation:

```
TMUE_LG2:19,449,820..19,460,071[rev] TTRE_chr2:12,163,850..12,172,620

TMUE_LG2:18,245,558..18,257,358[rev] TTRE_chr2:11,713,832..11,731,423 consider deleting the short tx

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
--->

NB: If you want to explore these data by yourself note that to keep the data
files small we only included genes and synteny links in these intervals:

*T. muris*: TMUE_LG2:17202488-20000820

*T. trichiura*: TTRE_chr2:10695245-13020708

