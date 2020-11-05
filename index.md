---
layout: page
title: Posts
tagline: blog
---

<div class="row">
  {% for post in site.posts limit:10 %}
  <div class="span3">

   <div style="float:right;">
    	{% if post.thumbnail %}
    	<img src="{{ post.thumbnail }}" class="thumbnail" />
    	{% else %}
    	<img src="/img/logo.png" class="thumbnail" />
    	{% endif %}
    </div>

    <h2>
     <a href="{{ BASE_PATH }}{{ post.url }}" style="text-decoration: none;">{{ post.title }}  </a>
     <h4>{{post.date | date_to_long_string}}</h4> </h2>


    {% assign paragraphs = post.content | split:'</h2>' %}
    {{ paragraphs[1] | strip_html | truncatewords:80 }}
    <p>
    <a class="btn" href="{{ BASE_PATH }}{{ post.url }}">Read more...</a>
    </p>
   </div>
   <hr/>

  {% endfor %}

</div>  

<div>
<h1> All posts </h1>
<ul class="posts">
  {% for post in site.posts %}
    <li><span>{{ post.date | date_to_string }}</span> &raquo; <a href="{{ BASE_PATH }}{{ post.url }}">{{ post.title }}</a></li>
  {% endfor %}
</ul>
</div>
