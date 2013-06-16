'use strict';

/***********************
 CLICKS / COMMON BROWSER EVENTS
 ***********************/

// note - requires jQuery
// future - to make more universal, see:
    // https://raw.github.com/cowboy/jquery-outside-events/v1.1/jquery.ba-outside-events.js
    // angular way of creating ng-directives
Application.Interface.directive('clickOutside', ['$document', '$parse', function($document, $parse) {
    return function(scope, element, attr, ctrl) {
        var handler = function(event) {
            //todo - rewrite without jQuery
            if (element.has(event.target).length == 0)
                scope.$apply(
                    $parse(attr['clickOutside'])(scope, {$event: event})
                );
        };

        $document.bind('click', handler);
        scope.$on('$destroy', function() {
            $document.unbind('click', handler);
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


Application.Interface.factory('keypressHelper', ['$parse', '$document', function keypress($parse, $document){
    var keysByCode = {
        8: 'backspace',
        9: 'tab',
        13: 'enter',
        27: 'esc',
        32: 'space',
        33: 'pageup',
        34: 'pagedown',
        35: 'end',
        36: 'home',
        37: 'left',
        38: 'up',
        39: 'right',
        40: 'down',
        45: 'insert',
        46: 'delete',
        188: 'comma',
        190: 'period'
    };

    var capitaliseFirstLetter = function (string) {
        return string.charAt(0).toUpperCase() + string.slice(1);
    };

    return function(mode, scope, elm, attrs) {
        var params, combinations = [];
        params = scope.$eval(attrs['ui'+capitaliseFirstLetter(mode)]);

        //CUSTOM
        if (elm == $document) {
            params = attrs;
        }

        // Prepare combinations for simple checking
        angular.forEach(params, function (v, k) {
            var combination, expression;
            expression = $parse(v);

            angular.forEach(k.split(' '), function(variation) {
                combination = {
                    expression: expression,
                    keys: {}
                };
                angular.forEach(variation.split('-'), function (value) {
                    combination.keys[value] = true;
                });
                combinations.push(combination);
            });
        });

        // Check only matching of pressed keys one of the conditions
        elm.bind(mode, function (event) {
            // No need to do that inside the cycle
            var metaPressed = !!(event.metaKey && !event.ctrlKey);
            var altPressed = !!event.altKey;
            var ctrlPressed = !!event.ctrlKey;
            var shiftPressed = !!event.shiftKey;
            var keyCode = event.keyCode;

            // normalize keycodes
            if (mode === 'keypress' && !shiftPressed && keyCode >= 97 && keyCode <= 122) {
                keyCode = keyCode - 32;
            }

            // Iterate over prepared combinations
            angular.forEach(combinations, function (combination) {

                var mainKeyPressed = combination.keys[keysByCode[event.keyCode]] || combination.keys[event.keyCode.toString()];

                var metaRequired = !!combination.keys.meta;
                var altRequired = !!combination.keys.alt;
                var ctrlRequired = !!combination.keys.ctrl;
                var shiftRequired = !!combination.keys.shift;

                if (
                    mainKeyPressed &&
                        ( metaRequired === metaPressed ) &&
                        ( altRequired === altPressed ) &&
                        ( ctrlRequired === ctrlPressed ) &&
                        ( shiftRequired === shiftPressed )
                    ) {
                    // Run the function
                    scope.$apply(function () {
                        combination.expression(scope, { '$event': event });
                    });
                }
            });
        });

        return [elm, mode, combinations];
    };
}]);

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


/***********************
 DIALOG / MODAL
 ***********************/

/* OLD VERSION --- NOT CLOSE TO COMPLETE
//todo - write - see angularStrap and AngularUI, bootstrapJS
//note - jQuery dependence
Application.Interface.directive('modal', ['Clotho', '$http', '$parse', '$templateCache', function(Clotho, $http, $parse, $templateCache) {
    return function(scope, element, attr, ctrl) {
        var background = $('<div class="modal-backdrop"></div>');

        //todo .... ummm?? async possible?
        var url = $parse(attr['modal'])(scope, {$event: event});
        scope.content = $http.get(url, {cache: $templateCache});

        element.bind('click', function(event) {
            scope.show(event);
        });

        scope.show = function(event) {
            background.appendTo('body').click(function(event) {
                scope.hide(event);
            });


        };

        scope.hide = function(event) {
            background.remove();
            event.preventDefault();
        }

    }
}]);*/

Application.Interface.provider("$dialog", function(){

    // The default options for all dialogs.
    var defaults = {
        backdrop: true,
        dialogClass: 'modal',
        backdropClass: 'modal-backdrop',
        transitionClass: 'fade',
        triggerClass: 'in',
        resolve:{},
        backdropFade: false,
        dialogFade:false,
        keyboard: true, // close with esc key
        backdropClick: true // only in conjunction with backdrop=true
        /* other options: template, templateUrl, controller */
    };

    var globalOptions = {};

    var activeBackdrops = {value : 0};

    // The `options({})` allows global configuration of all dialogs in the application.
    //
    //      var app = angular.module('App', ['ui.bootstrap.dialog'], function($dialogProvider){
    //        // don't close dialog when backdrop is clicked by default
    //        $dialogProvider.options({backdropClick: false});
    //      });
    this.options = function(value){
        globalOptions = value;
    };

    // Returns the actual `$dialog` service that is injected in controllers
    this.$get = ["$http", "$document", "$compile", "$rootScope", "$controller", "$templateCache", "$q", "$injector",
        function ($http, $document, $compile, $rootScope, $controller, $templateCache, $q, $injector) {

            var body = $document.find('body');

            function createElement(clazz) {
                var el = angular.element("<div>");
                el.addClass(clazz);
                return el;
            }

            // The `Dialog` class represents a modal dialog. The dialog class can be invoked by providing an options object
            // containing at lest template or templateUrl and controller:
            //
            //     var d = new Dialog({templateUrl: 'foo.html', controller: 'BarController'});
            //
            // Dialogs can also be created using templateUrl and controller as distinct arguments:
            //
            //     var d = new Dialog('path/to/dialog.html', MyDialogController);
            function Dialog(opts) {

                var self = this, options = this.options = angular.extend({}, defaults, globalOptions, opts);
                this._open = false;

                this.backdropEl = createElement(options.backdropClass);
                if(options.backdropFade){
                    this.backdropEl.addClass(options.transitionClass);
                    this.backdropEl.removeClass(options.triggerClass);
                }

                this.modalEl = createElement(options.dialogClass);
                if(options.dialogFade){
                    this.modalEl.addClass(options.transitionClass);
                    this.modalEl.removeClass(options.triggerClass);
                }

                this.handledEscapeKey = function(e) {
                    if (e.which === 27) {
                        self.close();
                        e.preventDefault();
                        self.$scope.$apply();
                    }
                };

                this.handleBackDropClick = function(e) {
                    self.close();
                    e.preventDefault();
                    self.$scope.$apply();
                };

                this.handleLocationChange = function() {
                    self.close();
                };
            }

            // The `isOpen()` method returns wether the dialog is currently visible.
            Dialog.prototype.isOpen = function(){
                return this._open;
            };

            // The `open(templateUrl, controller)` method opens the dialog.
            // Use the `templateUrl` and `controller` arguments if specifying them at dialog creation time is not desired.
            Dialog.prototype.open = function(templateUrl, controller){
                var self = this, options = this.options;

                if(templateUrl){
                    options.templateUrl = templateUrl;
                }
                if(controller){
                    options.controller = controller;
                }

                if(!(options.template || options.templateUrl)) {
                    throw new Error('Dialog.open expected template or templateUrl, neither found. Use options or open method to specify them.');
                }

                this._loadResolves().then(function(locals) {
                    var $scope = locals.$scope = self.$scope = locals.$scope ? locals.$scope : $rootScope.$new();

                    self.modalEl.html(locals.$template);

                    if (self.options.controller) {
                        var ctrl = $controller(self.options.controller, locals);
                        self.modalEl.children().data('ngControllerController', ctrl);
                    }

                    $compile(self.modalEl)($scope);
                    self._addElementsToDom();

                    // trigger tranisitions
                    setTimeout(function(){
                        if(self.options.dialogFade){ self.modalEl.addClass(self.options.triggerClass); }
                        if(self.options.backdropFade){ self.backdropEl.addClass(self.options.triggerClass); }
                    });

                    self._bindEvents();
                });

                this.deferred = $q.defer();
                return this.deferred.promise;
            };

            // closes the dialog and resolves the promise returned by the `open` method with the specified result.
            Dialog.prototype.close = function(result){
                var self = this;
                var fadingElements = this._getFadingElements();

                if(fadingElements.length > 0){
                    for (var i = fadingElements.length - 1; i >= 0; i--) {
                        //OLD WAY
                        //$transition(fadingElements[i], removeTriggerClass).then(onCloseComplete);
                        //NEW WAY
                        removeTriggerClass(fadingElements[i]);
                    }
                    onCloseComplete();
                    return;
                }

                this._onCloseComplete(result);

                function removeTriggerClass(el){
                    el.removeClass(self.options.triggerClass);
                }

                function onCloseComplete(){
                    if(self._open){
                        self._onCloseComplete(result);
                    }
                }
            };

            Dialog.prototype._getFadingElements = function(){
                var elements = [];
                if(this.options.dialogFade){
                    elements.push(this.modalEl);
                }
                if(this.options.backdropFade){
                    elements.push(this.backdropEl);
                }

                return elements;
            };

            Dialog.prototype._bindEvents = function() {
                if(this.options.keyboard){ body.bind('keydown', this.handledEscapeKey); }
                if(this.options.backdrop && this.options.backdropClick){ this.backdropEl.bind('click', this.handleBackDropClick); }

                this.$scope.$on('$locationChangeSuccess', this.handleLocationChange);
            };

            Dialog.prototype._unbindEvents = function() {
                if(this.options.keyboard){ body.unbind('keydown', this.handledEscapeKey); }
                if(this.options.backdrop && this.options.backdropClick){ this.backdropEl.unbind('click', this.handleBackDropClick); }
            };

            Dialog.prototype._onCloseComplete = function(result) {
                this._removeElementsFromDom();
                this._unbindEvents();

                //CUSTOM
                this.$scope.Clotho.silence(this.$scope.$id);

                this.deferred.resolve(result);
            };

            Dialog.prototype._addElementsToDom = function(){
                body.append(this.modalEl);

                if(this.options.backdrop) {
                    if (activeBackdrops.value === 0) {
                        body.append(this.backdropEl);
                    }
                    activeBackdrops.value++;
                }

                this._open = true;
            };

            Dialog.prototype._removeElementsFromDom = function(){
                this.modalEl.remove();

                if(this.options.backdrop) {
                    activeBackdrops.value--;
                    if (activeBackdrops.value === 0) {
                        this.backdropEl.remove();
                    }
                }
                this._open = false;
            };

            // Loads all `options.resolve` members to be used as locals for the controller associated with the dialog.
            Dialog.prototype._loadResolves = function(){
                var values = [], keys = [], templatePromise, self = this;

                //note: CUSTOM
                if (this.options.dependencies) {
                    var depPromise = Application.mixin(this.options.dependencies);
                    keys.push('dependencies');
                    values.push(depPromise);
                }

                if (this.options.template) {
                    templatePromise = $q.when(this.options.template);
                } else if (this.options.templateUrl) {
                    templatePromise = $http.get(this.options.templateUrl, {cache:$templateCache})
                        .then(function(response) { return response.data; });
                }

                angular.forEach(this.options.resolve || [], function(value, key) {
                    keys.push(key);
                    values.push(angular.isString(value) ? $injector.get(value) : $injector.invoke(value));
                });

                keys.push('$template');
                values.push(templatePromise);

                return $q.all(values).then(function(values) {
                    var locals = {};
                    angular.forEach(values, function(value, index) {
                        locals[keys[index]] = value;
                    });
                    locals.dialog = self;
                    return locals;
                });
            };

            // The actual `$dialog` service that is injected in controllers.
            return {
                // Creates a new `Dialog` with the specified options.
                dialog: function(opts){
                    return new Dialog(opts);
                },
                // creates a new `Dialog` tied to the default message box template and controller.
                //
                // Arguments `title` and `message` are rendered in the modal header and body sections respectively.
                // The `buttons` array holds an object with the following members for each button to include in the
                // modal footer section:
                //
                // * `result`: the result to pass to the `close` method of the dialog when the button is clicked
                // * `label`: the label of the button
                // * `cssClass`: additional css class(es) to apply to the button for styling
                messageBox: function(title, message, buttons){
                    return new Dialog({
                        templateUrl: 'interface/dialogMessagebox.html',
                        controller: 'MessageBoxController',
                        resolve:
                            {model: function() {
                                return {
                                    title: title,
                                    message: message,
                                    buttons: buttons
                                };
                            }
                    }});
                },

                serverAlert: function(message) {
                    return new Dialog({
                        backdrop: true,
                        backdropFade: true,
                        keyboard: true,
                        backdropClick: true,
                        templateUrl: 'interface/dialogMessagebox.html',
                        controller: 'ServerAlertController',
                        resolve:
                            {model: function() {
                                return {
                                    title: "Server Message",
                                    message: message,
                                    buttons: [{result:'ok', label: 'OK', cssClass: 'btn-primary'}]
                                };
                            }
                    }});
                }

            };
        }];
});

Application.Interface.controller('MessageBoxController', ['$scope', 'dialog', 'model', function($scope, dialog, model){
    $scope.title = model.title;
    $scope.message = model.message;
    $scope.buttons = model.buttons;
    $scope.close = function(res){
        dialog.close(res);
    };
}]);

Application.Interface.controller('ServerAlertController', ['$scope', 'dialog', 'model', 'Clotho', function($scope, dialog, model, Clotho){
    $scope.title = model.title;
    $scope.message = model.message;
    $scope.buttons = model.buttons;
    $scope.close = function(res){
        dialog.close(res);
    };

    //todo - better handling
    Clotho.listen('serverAlert', function() {
        $scope.close('Another alert appeared');
        Clotho.say($scope.message);
    }, $scope.$id);

}]);

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