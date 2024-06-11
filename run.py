import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable
from anadama2.tracked import TrackedDirectory

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Analysis Template"     #Update the description as needed
    ) 

# Setting additional custom arguments for workflow - run.py
workflow.add_argument(
    name="lines", 
    desc="Number of lines to trim [default: 10]", 
    default="10")

workflow.add_argument(
    name="metadata", 
    desc="Metadata for performing analysis [default: input/metadata.tsv]", 
    default="input/metadata.tsv")


# Parsing the workflow arguments
args = workflow.parse_args()

#Loading the config setting
args.config = 'etc/config.ini'

workflow.add_task(
    "alphadiversity_species.R",                            #Command 
    depends=[TrackedExecutable("alphadiversity_species.R")],                                 #Tracking executable dependencies
    targets=args.output,                                                       #Output target directory
    )                                                      #Additional arguments 

workflow.add_task(
    "extendeddatafig2b.R",                            #Command 
    depends=[TrackedExecutable("alphadiversity_species.R")],                                 #Tracking executable dependencies
    targets=args.output,                                                       #Output target directory
    )  

workflow.add_task(
    "figure1b_heatmaps.R",                            #Command 
    depends=[TrackedExecutable("alphadiversity_species.R")],                                 #Tracking executable dependencies
    targets=args.output,                                                       #Output target directory
    )  




#Task5 Add the document to the workflow
pdf_report=os.path.join(os.getcwd(),args.output,"pdfReport.pdf")
workflow.add_document(
    templates="doc/template.py",
    depends= [args.output+"/data.tsv.notabs"],
    targets=pdf_report, 
    vars={
        "introduction_text": "Demo Title"
    })

# Run the workflow
workflow.go()