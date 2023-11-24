---
layout: default
title: Vignettes
permalink: /vignettes
---

# Data analysis vignettes


<ul>
  {% for vignette in site.posts %}
  <li><a href="{{ vignette.url }}" class="vignette-preview">{{ vignette.title }}</a></li>
  {% endfor %}
</ul>