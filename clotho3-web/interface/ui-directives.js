'use strict';

/***********************
 CLICKS / COMMON BROWSER EVENTS
 ***********************/

/**
 * Taken from Angular UI
 *
 * General-purpose Event binding. Bind any event not natively supported by Angular
 * Pass an object with keynames for events to ui-event
 * Allows $event object and $params object to be passed
 *
 * @example <input ui-event="{ focus : 'counter++', blur : 'someCallback()' }">
 * @example <input ui-event="{ myCustomEvent : 'myEventHandler($event, $params)'}">
 *
 * @param ui-event {string|object literal} The event to bind to as a string or a hash of events with their callbacks
 */
Application.Interface.directive('uiEvent', ['$parse',
    function ($parse) {
        return function ($scope, elm, attrs) {
            var events = $scope.$eval(attrs.uiEvent);
            angular.forEach(events, function (uiEvent, eventName) {
                var fn = $parse(uiEvent);
                elm.bind(eventName, function (evt) {
                    var params = Array.prototype.slice.call(arguments);
                    //Take out first paramater (event object);
                    params = params.splice(1);
                    fn($scope, {$event: evt, $params: params});
                    if (!$scope.$$phase) {
                        $scope.$apply();
                    }
                });
            });
        };
    }]);

// note - requires jQuery
// future - to make more universal, see:
// https://raw.github.com/cowboy/jquery-outside-events/v1.1/jquery.ba-outside-events.js
// angular way of creating ng-directives
// note - Requires manual activiation to avoid too many $digests()
//      - activate with $scope.$broadcast('$event:$active')
Application.Interface.directive('clickOutside', ['$document', '$parse', function($document, $parse) {
    return function(scope, element, attr, ctrl) {
        var handler = function(event) {
            //todo - rewrite without jQuery.has()
            if (element.has(event.target).length == 0)
                scope.$safeApply(
                    $parse(attr['clickOutside'])(scope, {$event: event})
                );
        };

        //kill on scope desctruction
        scope.$on('$destroy', function() {
            $document.unbind('click', handler);
        });

        //custom events
        scope.$on("$event:$active", function(event, id) {
            if (scope.$id == id) {
                event.preventDefault();
                $document.bind('click', handler);
            }
        });

        scope.$on("$event:$inactive", function(event, id) {
            if (scope.$id == id) {
                event.preventDefault();
                $document.unbind('click', handler);
            }
        });
    }
}]);

//todo - test -- prevent on mouseleave child element -- prevent propagation
Application.Interface.directive('myMouseleave', ['$parse', function($parse) {
    return function(scope, element, attr) {
        element.bind('mouseleave', function(event) {
            event.stopPropagation();
            $parse(attr['myMouseleave'])(scope, {$event: event});
        });
    };
}]);

//future - remove polyfill once available-- currently not present in angular (1.1.5)
Application.Interface.directive('ngFocus', ['$parse', function($parse) {
    return function(scope, element, attr) {
        var fn = $parse(attr['ngFocus']);
        element.bind('focus', function(event) {
            scope.$apply(function() {
                fn(scope, {$event:event});
            });
        });
    }
}]);

//future - remove polyfill once available-- currently not present in angular (1.1.5)
Application.Interface.directive('ngBlur', ['$parse', function($parse) {
    return function(scope, element, attr) {
        var fn = $parse(attr['ngBlur']);
        element.bind('blur', function(event) {
            scope.$apply(function() {
                fn(scope, {$event:event});
            });
        });
    }
}]);

/***********************
 BUTTONS
 ***********************/

//from Angular UI Bootstrap
Application.Interface.directive('btnRadio', function () {
    var activeClass = 'active';
    var toggleEvent = 'click';

    return {
        require:'ngModel',
        link:function (scope, element, attrs, ngModelCtrl) {

            var value = scope.$eval(attrs.btnRadio);

            //model -> UI
            scope.$watch(function () {
                return ngModelCtrl.$modelValue;
            }, function (modelValue) {
                if (angular.equals(modelValue, value)){
                    element.addClass(activeClass);
                } else {
                    element.removeClass(activeClass);
                }
            });

            //ui->model
            element.bind(toggleEvent, function () {
                if (!element.hasClass(activeClass)) {
                    scope.$apply(function () {
                        ngModelCtrl.$setViewValue(value);
                    });
                }
            });
        }
    };
});

//from Angular UI Bootstrap
Application.Interface.directive('btnCheckbox', function() {

    var activeClass = 'active',
        toggleEvent = 'click';

    return {
        require:'ngModel',
        link:function (scope, element, attrs, ngModelCtrl) {

            var trueValue = scope.$eval(attrs.btnCheckboxTrue);
            var falseValue = scope.$eval(attrs.btnCheckboxFalse);

            trueValue = angular.isDefined(trueValue) ? trueValue : true;
            falseValue = angular.isDefined(falseValue) ? falseValue : false;

            //model -> UI
            scope.$watch(function () {
                return ngModelCtrl.$modelValue;
            }, function (modelValue) {
                if (angular.equals(modelValue, trueValue)) {
                    element.addClass(activeClass);
                } else {
                    element.removeClass(activeClass);
                }
            });

            //ui->model
            element.bind(toggleEvent, function () {
                scope.$apply(function () {
                    ngModelCtrl.$setViewValue(element.hasClass(activeClass) ? falseValue : trueValue);
                });
            });
        }
    };
});


/***********************
 KEYSTROKES
 ***********************/

/**
 * Bind one or more handlers to particular keys or their combination
 * @param hash {mixed} keyBindings Can be an object or string where keybinding expression of keys or keys combinations and AngularJS Exspressions are set. Object syntax: "{ keys1: expression1 [, keys2: expression2 [ , ... ]]}". String syntax: ""expression1 on keys1 [ and expression2 on keys2 [ and ... ]]"". Expression is an AngularJS Expression, and key(s) are dash-separated combinations of keys and modifiers (one or many, if any. Order does not matter). Supported modifiers are 'ctrl', 'shift', 'alt' and key can be used either via its keyCode (13 for Return) or name. Named keys are 'backspace', 'tab', 'enter', 'esc', 'space', 'pageup', 'pagedown', 'end', 'home', 'left', 'up', 'right', 'down', 'insert', 'delete', 'period', 'comma'.
 * @example <input ui-keypress="{enter:'x = 1', 'ctrl-shift-space':'foo()', 'shift-13':'bar()'}" /> <input ui-keypress="foo = 2 on ctrl-13 and bar('hello') on shift-esc" />
 **/
Application.Interface.directive('uiKeydown', ['keypressHelper', function(keypressHelper){
    return {
        link: function (scope, elm, attrs) {
            keypressHelper('keydown', scope, elm, attrs);
        }
    };
}]);

Application.Interface.directive('uiKeypress', ['keypressHelper', function(keypressHelper){
    return {
        link: function (scope, elm, attrs) {
            keypressHelper('keypress', scope, elm, attrs);
        }
    };
}]);

Application.Interface.directive('uiKeyup', ['keypressHelper', function(keypressHelper){
    return {
        link: function (scope, elm, attrs) {
            keypressHelper('keyup', scope, elm, attrs);
        }
    };
}]);

Application.Interface.directive('restrictInput', ['$parse', function($parse) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, iElement, iAttrs, controller) {
            scope.$watch(iAttrs.ngModel, function(value) {
                if (!value) {
                    return;
                }
                $parse(iAttrs.ngModel).assign(scope, value.replace(new RegExp(iAttrs.restrict, 'g'), '').replace(/\s+/g, '-'));
            });
        }
    }
}]);

/***********************
 DIALOG / MODAL
 ***********************/

Application.Interface.directive('modal', ['$parse', '$dialog', function($parse, $dialog) {
    return {
        restrict: 'EA',
        terminal: true,
        link: function(scope, elm, attrs) {
            var opts = angular.extend({}, scope.$eval(attrs.options));
            var shownExpr = attrs.modal || attrs.show;
            var setClosed;

            // Create a dialog with the template as the contents of the directive
            // Add the current scope as the resolve in order to make the directive scope as a dialog controller scope
            opts = angular.extend(opts, {
                template: elm.html(),
                resolve: { $scope: function() { return scope; } }
            });
            var dialog = $dialog.dialog(opts);

            elm.remove();

            if (attrs.close) {
                setClosed = function() {
                    $parse(attrs.close)(scope);
                };
            } else {
                setClosed = function() {
                    if (angular.isFunction($parse(shownExpr).assign)) {
                        $parse(shownExpr).assign(scope, false);
                    }
                };
            }

            scope.$watch(shownExpr, function(isShown, oldShown) {
                if (isShown) {
                    dialog.open().then(function(){
                        setClosed();
                    });
                } else {
                    //Make sure it is not opened
                    if (dialog.isOpen()){
                        dialog.close();
                    }
                }
            });
        }
    };
}]);


/******************
 TOOLTIP
 ******************/

/**
 * The following features are still outstanding:
 * animation as a function
 * placement as a function
 * inside
 * support for more triggers than just mouse enter/leave
 * html tooltips
 * and selector delegation.
 */

/**
 * The $tooltip service creates tooltip- and popover-like directives as well as
 * houses global options for them.
 */
Application.Interface.provider( '$tooltip', function () {
    // The default options tooltip and popover.
    var defaultOptions = {
        placement: 'top',
        animation: true,
        popupDelay: 0
    };

    // Default hide triggers for each show trigger
    var triggerMap = {
        'mouseenter': 'mouseleave',
        'click': 'click',
        'focus': 'blur'
    };

    // The options specified to the provider globally.
    var globalOptions = {};

    /**
     * `options({})` allows global configuration of all tooltips in the
     * application.
     *
     *   var app = angular.module( 'App', ['ui.bootstrap.tooltip'], function( $tooltipProvider ) {
   *     // place tooltips left instead of top by default
   *     $tooltipProvider.options( { placement: 'left' } );
   *   });
     */
    this.options = function( value ) {
        angular.extend( globalOptions, value );
    };

    /**
     * This allows you to extend the set of trigger mappings available. E.g.:
     *
     *   $tooltipProvider.setTriggers( 'openTrigger': 'closeTrigger' );
     */
    this.setTriggers = function setTriggers ( triggers ) {
        angular.extend( triggerMap, triggers );
    };

    /**
     * This is a helper function for translating camel-case to snake-case.
     */
    function snake_case(name){
        var regexp = /[A-Z]/g;
        var separator = '-';
        return name.replace(regexp, function(letter, pos) {
            return (pos ? separator : '') + letter.toLowerCase();
        });
    }

    /**
     * Returns the actual instance of the $tooltip service.
     */
    this.$get = [ '$window', '$compile', '$timeout', '$parse', '$document', '$position', '$interpolate', function ( $window, $compile, $timeout, $parse, $document, $position, $interpolate ) {
        return function $tooltip ( type, prefix, defaultTriggerShow ) {
            var options = angular.extend( {}, defaultOptions, globalOptions );

            /**
             * Returns an object of show and hide triggers.
             *
             * If a trigger is supplied,
             * it is used to show the tooltip; otherwise, it will use the `trigger`
             * option passed to the `$tooltipProvider.options` method; else it will
             * default to the trigger supplied to this directive factory.
             *
             * The hide trigger is based on the show trigger. If the `trigger` option
             * was passed to the `$tooltipProvider.options` method, it will use the
             * mapped trigger from `triggerMap` or the passed trigger if the map is
             * undefined; otherwise, it uses the `triggerMap` value of the show
             * trigger; else it will just use the show trigger.
             */
            function setTriggers ( trigger ) {
                var show, hide;

                show = trigger || options.trigger || defaultTriggerShow;
                if ( angular.isDefined ( options.trigger ) ) {
                    hide = triggerMap[options.trigger] || show;
                } else {
                    hide = triggerMap[show] || show;
                }

                return {
                    show: show,
                    hide: hide
                };
            }

            var directiveName = snake_case( type );
            var triggers = setTriggers( undefined );

            var startSym = $interpolate.startSymbol();
            var endSym = $interpolate.endSymbol();
            var template =
                '<'+ directiveName +'-popup '+
                    'title="'+startSym+'tt_title'+endSym+'" '+
                    'content="'+startSym+'tt_content'+endSym+'" '+
                    'placement="'+startSym+'tt_placement'+endSym+'" '+
                    'animation="tt_animation()" '+
                    'is-open="tt_isOpen"'+
                    '>'+
                    '</'+ directiveName +'-popup>';

            return {
                restrict: 'EA',
                scope: true,
                link: function link ( scope, element, attrs ) {
                    var tooltip = $compile( template )( scope );
                    var transitionTimeout;
                    var popupTimeout;
                    var $body;
                    var appendToBody = angular.isDefined( options.appendToBody ) ? options.appendToBody : false;

                    // By default, the tooltip is not open.
                    // TODO add ability to start tooltip opened
                    scope.tt_isOpen = false;

                    function toggleTooltipBind () {
                        if ( ! scope.tt_isOpen ) {
                            showTooltipBind();
                        } else {
                            hideTooltipBind();
                        }
                    }

                    // Show the tooltip with delay if specified, otherwise show it immediately
                    function showTooltipBind() {
                        if ( scope.tt_popupDelay ) {
                            popupTimeout = $timeout( show, scope.tt_popupDelay );
                        } else {
                            scope.$apply( show );
                        }
                    }

                    function hideTooltipBind () {
                        scope.$apply(function () {
                            hide();
                        });
                    }

                    // Show the tooltip popup element.
                    function show() {
                        var position,
                            ttWidth,
                            ttHeight,
                            ttPosition;

                        // Don't show empty tooltips.
                        if ( ! scope.tt_content ) {
                            return;
                        }

                        // If there is a pending remove transition, we must cancel it, lest the
                        // tooltip be mysteriously removed.
                        if ( transitionTimeout ) {
                            $timeout.cancel( transitionTimeout );
                        }

                        // Set the initial positioning.
                        tooltip.css({ top: 0, left: 0, display: 'block' });

                        // Now we add it to the DOM because need some info about it. But it's not
                        // visible yet anyway.
                        if ( appendToBody ) {
                            $body = $body || $document.find( 'body' );
                            $body.append( tooltip );
                        } else {
                            element.after( tooltip );
                        }

                        // Get the position of the directive element.
                        position = appendToBody ? $position.offset( element ) : $position.position( element );

                        // Get the height and width of the tooltip so we can center it.
                        ttWidth = tooltip.prop( 'offsetWidth' );
                        ttHeight = tooltip.prop( 'offsetHeight' );

                        // Calculate the tooltip's top and left coordinates to center it with
                        // this directive.
                        switch ( scope.tt_placement ) {
                            case 'mouse':
                                var mousePos = $position.mouse();
                                ttPosition = {
                                    top: mousePos.y,
                                    left: mousePos.x
                                };
                                break;
                            case 'right':
                                ttPosition = {
                                    top: position.top + position.height / 2 - ttHeight / 2,
                                    left: position.left + position.width
                                };
                                break;
                            case 'bottom':
                                ttPosition = {
                                    top: position.top + position.height,
                                    left: position.left + position.width / 2 - ttWidth / 2
                                };
                                break;
                            case 'left':
                                ttPosition = {
                                    top: position.top + position.height / 2 - ttHeight / 2,
                                    left: position.left - ttWidth
                                };
                                break;
                            default:
                                ttPosition = {
                                    top: position.top - ttHeight,
                                    left: position.left + position.width / 2 - ttWidth / 2
                                };
                                break;
                        }

                        ttPosition.top += 'px';
                        ttPosition.left += 'px';

                        // Now set the calculated positioning.
                        tooltip.css( ttPosition );

                        // And show the tooltip.
                        scope.tt_isOpen = true;
                    }

                    // Hide the tooltip popup element.
                    function hide() {
                        // First things first: we don't show it anymore.
                        scope.tt_isOpen = false;

                        //if tooltip is going to be shown after delay, we must cancel this
                        $timeout.cancel( popupTimeout );

                        // And now we remove it from the DOM. However, if we have animation, we
                        // need to wait for it to expire beforehand.
                        // FIXME: this is a placeholder for a port of the transitions library.
                        if ( angular.isDefined( scope.tt_animation ) && scope.tt_animation() ) {
                            transitionTimeout = $timeout( function () { tooltip.remove(); }, 500 );
                        } else {
                            tooltip.remove();
                        }
                    }

                    /**
                     * Observe the relevant attributes.
                     */
                    attrs.$observe( type, function ( val ) {
                        scope.tt_content = val;
                    });

                    attrs.$observe( prefix+'Title', function ( val ) {
                        scope.tt_title = val;
                    });

                    attrs.$observe( prefix+'Placement', function ( val ) {
                        scope.tt_placement = angular.isDefined( val ) ? val : options.placement;
                    });

                    attrs.$observe( prefix+'Animation', function ( val ) {
                        scope.tt_animation = angular.isDefined( val ) ? $parse( val ) : function(){ return options.animation; };
                    });

                    attrs.$observe( prefix+'PopupDelay', function ( val ) {
                        var delay = parseInt( val, 10 );
                        scope.tt_popupDelay = ! isNaN(delay) ? delay : options.popupDelay;
                    });

                    attrs.$observe( prefix+'Trigger', function ( val ) {
                        element.unbind( triggers.show );
                        element.unbind( triggers.hide );

                        triggers = setTriggers( val );

                        if ( triggers.show === triggers.hide ) {
                            element.bind( triggers.show, toggleTooltipBind );
                        } else {
                            element.bind( triggers.show, showTooltipBind );
                            element.bind( triggers.hide, hideTooltipBind );
                        }
                    });

                    attrs.$observe( prefix+'AppendToBody', function ( val ) {
                        appendToBody = angular.isDefined( val ) ? $parse( val )( scope ) : appendToBody;
                    });

                    // if a tooltip is attached to <body> we need to remove it on
                    // location change as its parent scope will probably not be destroyed
                    // by the change.
                    if ( appendToBody ) {
                        scope.$on('$locationChangeSuccess', function closeTooltipOnLocationChangeSuccess () {
                            if ( scope.tt_isOpen ) {
                                hide();
                            }
                        });
                    }

                    // Make sure tooltip is destroyed and removed.
                    scope.$on('$destroy', function onDestroyTooltip() {
                        if ( scope.tt_isOpen ) {
                            hide();
                        } else {
                            tooltip.remove();
                        }
                    });
                }
            };
        };
    }];
});

Application.Interface.directive( 'tooltipPopup', function () {
    return {
        restrict: 'E',
        replace: true,
        scope: { content: '@', placement: '@', animation: '&', isOpen: '&' },
        templateUrl: 'interface/templates/tooltip-popup.html'
    };
});

Application.Interface.directive( 'tooltip', [ '$tooltip', function ( $tooltip ) {
    return $tooltip( 'tooltip', 'tooltip', 'mouseenter' );
}]);

Application.Interface.directive( 'tooltipHtmlUnsafePopup', function () {
    return {
        restrict: 'E',
        replace: true,
        scope: { content: '@', placement: '@', animation: '&', isOpen: '&' },
        templateUrl: 'interface/templates/tooltip-html-unsafe-popup.html'
    };
});

Application.Interface.directive( 'tooltipHtmlUnsafe', [ '$tooltip', function ( $tooltip ) {
    return $tooltip( 'tooltipHtmlUnsafe', 'tooltip', 'mouseenter' );
}]);