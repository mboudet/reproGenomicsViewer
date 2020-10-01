import os
import json
from datetime import datetime

from dal import autocomplete
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic
from django.utils import timezone
from django.views import View
from django.shortcuts import redirect
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.http import JsonResponse
from django.contrib import messages
from django.template.loader import render_to_string
from django.db.models import Q
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.generic import CreateView
from django.core import serializers

from guardian.shortcuts import get_perms

import pandas as pd
import numpy as np
import json
import statistics

import uuid
import shutil
from .models import ExpressionStudy, ExpressionData, Gene, Database
from .forms import *
from .graphs import getClasses, get_graph_data_full, get_graph_data_genes, getValues, getValuesExpression, get_density_graph_data_full, get_density_graph_gene_data_full, get_violin_graph_gene_data_full

class GeneAutocomplete(autocomplete.Select2QuerySetView):

    def get_result_value(self, result):
        return result.id

    def get_result_label(self, result):
        return result.symbol

    def get_queryset(self):
        query = self.q
        qs = Gene.objects.all()
        if query:
            qs = qs.filter(Q(symbol__icontains=query) | Q(synonyms__icontains=query)| Q(gene_id__icontains=query))
        return qs

def get_gene(request, gene_id):

    gene = get_object_or_404(Gene, id=gene_id)
    data = {'gene_id' : gene.gene_id, "symbol": gene.symbol, "homolog_id": gene.homolog_id, "ensembl_id": gene.ensemble_id}
    return JsonResponse(data)

def get_stud_db(request):
    db_ids  = request.GET.getlist('db_ids[]')    
    print (db_ids)
    db_ids = list(map(int, db_ids))
    print (db_ids)
    studies_list = []
    for db in db_ids :
        print(db)
        database =get_object_or_404(Database, id=db)
        print(database)
        studies = ExpressionStudy.objects.filter(database=database)
        studies_list.extend(studies)

    print(studies_list)

    table = render_to_string('studies/partial_study_table.html', {'studies': studies_list}, request)
    data = {'table_list' : table}

    return JsonResponse(data)

def index(request):
    columns = [
            "article",
            "pmid",
            "ome",
            "technology",
            "species",
            "experimental_design",
            "topics",
            "tissues",
            "sex",
            "dev_stage",
            "age",
            "antibody",
            "mutant",
            "cell_sorted",
            "keywords",
    ]

    all_studies = [study for study in ExpressionStudy.objects.exclude(data=None) if check_view_permissions(request.user, study)]
    studies = paginate(all_studies)
    form = ExpressionStudyFilterForm(studies=all_studies)
    table = render_to_string('studies/partial_study_table.html', {'studies': studies}, request)
    pagination = render_to_string('studies/partial_study_pagination.html', {'table': studies}, request)
    context = {'form': form, 'columns': columns, 'table': table, 'pagination': pagination}
    return render(request, 'studies/scatter_plot.html', context)

def document_select(request):

    if not "id" in request.GET:
        return redirect(reverse("studies:index"))

    id_list = request.GET.getlist("id")
    # Just in case
    if not all(x.isdigit() for x in id_list):
        return redirect(reverse("studies:index"))

    studies = ExpressionStudy.objects.filter(id__in=id_list)
    if studies.count() == 0:
        return redirect(reverse("studies:index"))

    table = render_to_string('studies/document_select.html', {'studies': studies}, request)
    data = {'table' : table}
    
    return JsonResponse(data)
    #return render(request, 'studies/document_select.html', {'studies': studies})

def show_graph(request):

    if not "document_id" in request.GET and not "study_id" in request.GET:
        return redirect(reverse("studies:index"))

    document_id = request.GET["document_id"]
    study_id = request.GET["study_id"]
    # Just in case
    if not document_id.isdigit() or not study_id.isdigit():
        return redirect(reverse("studies:index"))


    ##########################
    # File Statistique only for RGV files
    ##########################
    data = get_object_or_404(ExpressionData, id=document_id)

    data_stat = {}

        
    study = get_object_or_404(ExpressionStudy, id=study_id)
    form = GeneFilterForm()
    classes = getClasses(data)
    context = {'study': study, 'document': data, 'classes': classes, 'form': form}
    return render(request, 'studies/graph.html', context)

def get_graph_data(request):
    
    display_mode = "scatter"

    if "mode" in request.GET:
        display_mode = request.GET["mode"]
    

    if not "document_id" in request.GET:
        return redirect(reverse("studies:index"))

    document_id = request.GET["document_id"]

    if not document_id.isdigit():
        return redirect(reverse("studies:index"))

    data = get_object_or_404(ExpressionData, id=document_id)
    

    selected_class = request.GET.get('selected_class', None)

    if "gene_id" in request.GET:

        exp_list = []
        if "|" in request.GET['gene_id']:
            list_gene = request.GET['gene_id'].split("|")
            for g_ in list_gene :
                gene = get_object_or_404(Gene, id=g_)
                exp_list.append(gene)
        else :
            gene = get_object_or_404(Gene, id=request.GET['gene_id'])
            exp_list.append(gene)
        if display_mode =="scatter" :
            data = get_graph_data_genes(data,exp_list, selected_class)
        if display_mode =="density" :
            data = get_density_graph_gene_data_full(data,exp_list, selected_class)
        if display_mode =="violin" :
            data = get_violin_graph_gene_data_full(data,exp_list, selected_class)
    else:
        if display_mode =="scatter" :
            data = get_graph_data_full(data, selected_class)
        if display_mode =="density" :
            data = get_density_graph_data_full(data, selected_class)
    return JsonResponse(data)

def get_group_info(request):

    group = request.GET.get('group',None)
    sample = request.GET.get('sample',None)
    document_id = request.GET.get('document',None)
    
    data = get_object_or_404(ExpressionData, id=document_id)

    selected_class = request.GET.get('selected_class', None)

    expression_values_group = getValuesExpression(data, selected_class, group)


    table = render_to_string('studies/partial_group_info.html', {'genes_list': expression_values_group}, request)
    data = {'list' : expression_values_group, 'group':group}



    return JsonResponse(data)


def render_table(request):

    data = {}
    studies = ExpressionStudy.objects.exclude(data=None)
    kwargs = {}
    for key, value in request.GET.items():
        if value:
            if key == "article":
                kwargs[key + "__icontains"] = value
            elif key == "technology" or key == "species":
                kwargs["data__" + key] = value
            elif key == "page":
                continue
            else:
                kwargs[key + "__contains"] = [value]

    
    studies = paginate([study for study in studies.filter(**kwargs).distinct() if check_view_permissions(request.user, study)], request.GET.get('page'))
    # Filter here
    table = render_to_string('studies/partial_study_table.html', {'studies': studies}, request)
    pagination = render_to_string('studies/partial_study_pagination.html', {'table': studies}, request)
    data['table'] = table
    data['pagination'] = pagination
    return JsonResponse(data)

def autocomplete_genes(request,taxonid):

    if request.is_ajax():
        query = request.GET.get('term','')
        qs = Gene.objects.all()
        qs = qs.filter(Q(tax_id__exact=int(taxonid)))
        qs = qs.filter(Q(symbol__icontains=query) | Q(synonyms__icontains=query)| Q(gene_id__icontains=query) & Q(tax_id__exact=int(taxonid)))
        results = []
        for gene in qs :
            results.append({'label' : gene.symbol, 'value':gene.symbol+" ("+str(gene.id)+")"})
        data = json.dumps(results[:10])
    else:
        data="fail"
    mimetype = 'application/json'
    return HttpResponse(data, mimetype)

def paginate(values, query=None, count=5, is_ES=False):

    paginator = Paginator(values, count)

    try:
        val = paginator.page(query)
    except PageNotAnInteger:
        val = paginator.page(1)
    except EmptyPage:
        val = paginator.page(paginator.num_pages)

    return val

def check_view_permissions(user, study, strict=False):
    has_access = False
    if study.status == "PUBLIC" and not strict:
        has_access = True
    elif user.is_superuser:
        has_access = True
    elif user.is_authenticated and 'view_expressionstudy' in get_perms(user, study):
        has_access = True

    return has_access
