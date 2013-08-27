# Clotho 3.0 Sharable Editor - AngularJS

## Overview

Clotho 3.0 Front-end Playground: Functionalities for:
- Sharable Editor
- Trails Browser
- Chat (socket implementation)

## To Run

- set current directory to this one, then run `node scripts/web-server-socket.js` and go to localhost:8080/app/index.html
- note: socket.io must be installed locally. See http://socket.io/#how-to-use

## Prerequisites

### Git

### Node.js
- web-server currently runs on node.

### socket.io
- Must be installed for Socket communication. See: http://socket.io/#how-to-use

### Java
- http://www.java.com




## References

### From Angular

- DOCS
- API

#### Notes
- https://github.com/angular/angular.js/wiki/The-Nuances-of-Scope-Prototypal-Inheritance

### Personal notes

#### Compatibility
- avoid directive element names (type E) for IE compatibility, or use document.createElement (i.e. for non-native elements like <widget></widget>)
    - see AngularUI ie-shiv


## TODOs, Future references, etc.

### Important Fixes + Changes



### Current Tasks

- anderson lab site

- finish schema editor -- need stephanie response

- search bar rewrite
    - autocomplete positioning using $position

- login window

- API
    - Socket events in own layer (outside PubSub)
    - verify collector updates on set() and destroy() and get()
        - collect queries
    - add-in $dialog for alert() and edit() but don't make dependency
        - remove from anderson_lab dependencies


- Trails
    - start writing demo trail
    - Trail schema writeup / Trail
    - Trail service to persist variables -- Trails.persist()
    - Trail editor -- req. schema finalized
    - exercise template (intro, tempalate / video, question)
    - prettyprint directive
    - start approaching construction files


- Random Fixes
    - search module requests [object Object] -- find where and remove

- Mac Housekeeping
    - folder actions to push changes to models folder to testData folder
    - tie github pulls to original repositories and keep in sync
    - package clotho api stuff (foundation folder and application.js)
        - grunt / yeoman --- javascript project builder


- DNA Functions
    - decide data structure
    - DNA
        - required upstream
            - determine ORF
            - species specific
                - codon optimize
                - silent sites
        - parse: gb, fasta, embl
        - BLAST breakout
    - Digest
        - import NEB enzymes, use same format
            - call NEB @ 1-800-632-7799
            - cuts: ^   |   _
            - broaden methylation categories: dam, dcm, cpg
        - dep on DNA service new functions
            - swapout Digest site (need to know ORF)
        - import strider format (REBASE) / NEB format
    - PCR
        - other pcr types


- editor
    - replace internal HTML, not whole thing (closeable directive)
    - use "ID" not "UUID"
    - dependencies field should query other functions
    - use accordion for tests div
    - map BSON -> HTML5 data types
    - query for test args
        - show name, not uuid
        - Popover with whole JSON shown
    - use form controller for validation
        - currently requires nested forms: https://github.com/angular/angular.js/issues/1404
    - new fields
        - multiple select
        - select
        - radios
        - file drops?
    - add novalidate via directive
        - need to test novalidate="true"
        - custom validation using ngForm or ngModel.$setValidity

    - function editor
        - integrate codemirror, tie to language selection


- UI Directives
    - popover with HTML (see angularStrap, or wait for bootstrap UI)
        - https://github.com/angular-ui/bootstrap/issues/220
    - rewrite contenteditable
        - make play nicely (e.g. with digest-highlight)
    - dropdown in Angular --  see angular-strap
    - popover with HTML -- test -- see angular-strap
        - hide on mouseleave when mouseenter trigger
    - angular Strap: https://raw.github.com/mgcrea/angular-strap/v0.7.5/dist/angular-strap.js
    - DRAGGABLE
        - https://github.com/fisshy/Angular-drag-drop
        - need to handle drop + other events


- FUTURE : move to angular 1.2.0


- plasmid editors to check out
    - benchling
    - sparkDNA



- caret movement directives : https://github.com/DrPheltRight/jquery-caret
    - see also contenteditable directive - have one working
    - get current position and save as option



- bootstrapping apps without ng-view (just use ng-include or templates)
    - use ng 1.2.0 - routes not included by default



- create wrappers for directives, filters -> expose DNA functions accordingly



- break up UI directives, filters, etc. into modules



- draggable elements (within angular -- also sortable)
    - https://github.com/fisshy/Angular-drag-drop/blob/master/javascript/dragandrop.js
    - http://jsfiddle.net/ADukg/2516/



- (advantage?) plasmid editor -- add locations parser to ngModel pipeline, store locations inside of ngModel (on the object), push them back in using a formatter



- $keypress.suspend() e.g. in trail mode when using search bar -- need to be namespaced


- move searchbar to typeahead


- trails functionality
    - youtube events - NOT WORKING
        - pass as param to function
    - templates
        - types
            - sequence
            - assertion / reason
            - hot spot (choose the wrong thing) /// paragraph with multi TFs
            - multipic - add active class
            - drag-drop - make draggable
            - ranking - tie to model
    - UI directives
        - draggable
            - see UI-sortable
            - drag-to for drag-drop
    - answering quiz colors + disables side-nav


- Search Bar
    - reset autocomplete detail when changing
    - animate showing + hiding -- use $animate
    - autocomplete
        - use cache for results?
        - $timeout sending to 50ms?
    - mouseenter and leave events -- maybe size list items differently
        - undetail propagate to children? why?
        - see: http://plnkr.co/edit/QSgRu4Sr8TSqR86eAfvq?p=preview
    - future: rewrite
        - easy to transclude, with customizable / default options



- favorite() API method - and unfavorite() and dislike()


- clickOutside handlers - better logic esp. for ID


- Dynamic Adds
    - Reference
        - watch : https://github.com/mhevery/angular.js/blob/dte/src/auto/injector.js
        - controllers : https://github.com/matys84pl/angularjs-requirejs-lazy-controllers




UI STUFF
- services to write / borrow
    - dropdown (bootstrap)
        - use for searchbar help
    - progressBar (bootstrap)
    - scrollTo (probably easy to write) / scrollfix (see angular ui utils)
- ANIMATION
    - see CSS framework https://github.com/daneden/animate.css


- Clotho.show()
    - store model and template name somewhere, access in dynamic_template for show()
        - after download JS: callback or something after controller registered?
        - Register as controller object etc. rather than global function
            - EXAMPLE SEE http://stackoverflow.com/questions/15250644/angularjs-loading-a-controller-dynamically
        - Caching templates
            - remove and reset if exists
                - http://deansofer.com/posts/view/14/AngularJs-Tips-and-Tricks-UPDATED#refresh-partials


    - show() in modal box as option --- make this alert?
        - simple example: http://jsfiddle.net/jjzabkar/MWAek/1/
        - see angular UI / UI-bootstrap modal directive
            - http://angular-ui.github.io/bootstrap/




- TESTING
    - write unit tests, etc.




- Clotho API
    - once REST works, write synchronous Clotho.get()
    - combine watch and watch2 - check for function vs. obj and field
        - try using angular.extend()



- remove listeners
    - see $scope.$watch:
        - deregistration returned as function to invoke later
        - http://www.bennadel.com/blog/2480-Unbinding-watch-Listeners-In-AngularJS.htm
    - reg app ID with each app
    - make universal --> $scope.$on('$destroy', Clotho.silence($scope.$id))
            -i.e. (?) extend $rootScope.$destroy() if possible



- Refine PubSub / communication
    - REMOVE SOCKET LISTENERS (only listen to PubSub)
    - refine pubsub listeners object
    - PubSub clear needs to clean 'references' object
    - references
        - http://closure-library.googlecode.com/svn/docs/closure_goog_pubsub_pubsub.js.source.html
        - http://www.gridlinked.info/angularJS/modules/MessagingServices.js (requires jQuery)



- Routing
    - app with no routes (i.e. for widgets) ???
    - test asynchronously adding routes via controller
        - reference
           - http://stackoverflow.com/questions/13153121/how-to-defer-routes-definition-in-angular-js
           - https://groups.google.com/forum/#!msg/angular/mrcy_2BZavQ/Mqte8AvEh0QJ
           - http://stackoverflow.com/questions/13681116/angularjs-dynamic-routing



- deep-linking for search bar and editor
    - http://stackoverflow.com/questions/14974271/can-you-change-a-path-without-reloading-the-controller-in-angularjs





---------- maybes --------------


- (?) write ClothoProvider - include in module config
    - look for other advantages / utility



- decorators? i.e. for $http


------ experiments ------

- filters that perform an action / broadcast events / etc.



#### Research

- re-read: http://deansofer.com/posts/view/14/AngularJs-Tips-and-Tricks-UPDATED#refresh-partials
- currying? angular.bind? use cases?

### Future tasks

- D3.js to make simple plasmids

### TODOs by category


#### General UI

- alerts via Bootstrap (and modals?)

#### General Backend


#### Forward looking


### At Deployment

#### General

- minify and concatenate JS
- non-cacheable template
- cache with versions: js, css, views, code
- sanitize (escape) HTML
- IE <7 compatibility
    - document.createElement(<custom tag names>)
    - JSON polyfill

#### AJAX Crawlable

- https://developers.google.com/webmasters/ajax-crawling/docs/specification

## Things / Apps to check out

#### General Examples

- https://github.com/angular/angular.js/wiki/JsFiddle-Examples
- http://builtwith.angularjs.org/

#### ActiveMQ & node.js

- http://stackoverflow.com/questions/4700935/what-are-good-message-queue-options-for-nodejs

#### PubNub

- http://cdn.pubnub.com/pubnub-3.1.js
- http://jsfiddle.net/bv5Kq/12/

#### Authentication

- http://www.espeo.pl/2012/02/26/authentication-in-angularjs-application

#### Widgets

- http://jsfiddle.net/simpulton/VJ94U/
- http://onehungrymind.com/angularjs-sticky-notes-pt-1-architecture/
- Simple Draggable Directive: http://docs.angularjs.org/guide/compiler
- Draggable things in Angular (with shared controller w/o require): http://www.bennadel.com/blog/2446-Using-Controllers-In-Directives-In-AngularJS.htm

#### Manual Bootstrapping / Multiple Apps per page

- http://docs.angularjs.org/api/angular.bootstrap
- https://groups.google.com/forum/#!msg/angular/lhbrIG5aBX4/4hYnzq2eGZwJ