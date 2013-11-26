# Clotho 3.0

-- add in WebSocket as part of config block so ready before anything runs

## Todo

LOW PRIORITY

- UI services
    - nice dropdown for filtering
          http://jsfiddle.net/TahmidTanzim/N9Vqk/


- upgrade to $modal from angular ui


- UI Directives
    - popover with HTML (see angularStrap, or wait for bootstrap UI)
        - https://github.com/angular-ui/bootstrap/issues/220
        - hide on mouseleave when mouseenter trigger
        - angular Strap: https://raw.github.com/mgcrea/angular-strap/v0.7.5/dist/angular-strap.js
    - rewrite contenteditable
        - make play nicely (e.g. with digest-highlight)


- Clotho API
    - once REST works, write synchronous Clotho.get()

------ experiments ------

- (cool idea // advantage?) plasmid editor -- add locations parser to ngModel pipeline, store locations inside of ngModel (on the object), push them back in using a formatter

- filters that perform an action / broadcast events / etc.
