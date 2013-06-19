'use strict';

/**
 * taken from Angular UI Bootstrap
 * note: only required if jQuery not present
 *
 * A set of utility methods that can be use to retrieve position of DOM elements.
 * It is meant to be used where we need to absolute-position DOM elements in
 * relation to other, existing elements (this is the case for tooltips, popovers,
 * typeahead suggestions etc.).
 */

Application.Interface.factory('$position', ['$document', '$window', function ($document, $window) {

    function getStyle(el, cssprop) {
        if (el.currentStyle) { //IE
            return el.currentStyle[cssprop];
        } else if ($window.getComputedStyle) {
            return $window.getComputedStyle(el)[cssprop];
        }
        // finally try and get inline style
        return el.style[cssprop];
    }

    /**
     * Checks if a given element is statically positioned
     * @param element - raw DOM element
     */
    function isStaticPositioned(element) {
        return (getStyle(element, "position") || 'static' ) === 'static';
    }

    /**
     * returns the closest, non-statically positioned parentOffset of a given element
     * @param element
     */
    var parentOffsetEl = function (element) {
        var docDomEl = $document[0];
        var offsetParent = element.offsetParent || docDomEl;
        while (offsetParent && offsetParent !== docDomEl && isStaticPositioned(offsetParent) ) {
            offsetParent = offsetParent.offsetParent;
        }
        return offsetParent || docDomEl;
    };

    return {
        /**
         * Provides read-only equivalent of jQuery's position function:
         * http://api.jquery.com/position/
         */
        position: function (element) {
            var elBCR = this.offset(element);
            var offsetParentBCR = { top: 0, left: 0 };
            var offsetParentEl = parentOffsetEl(element[0]);
            if (offsetParentEl != $document[0]) {
                offsetParentBCR = this.offset(angular.element(offsetParentEl));
                offsetParentBCR.top += offsetParentEl.clientTop;
                offsetParentBCR.left += offsetParentEl.clientLeft;
            }

            return {
                width: element.prop('offsetWidth'),
                height: element.prop('offsetHeight'),
                top: elBCR.top - offsetParentBCR.top,
                left: elBCR.left - offsetParentBCR.left
            };
        },

        /**
         * Provides read-only equivalent of jQuery's offset function:
         * http://api.jquery.com/offset/
         */
        offset: function (element) {
            var boundingClientRect = element[0].getBoundingClientRect();
            return {
                width: element.prop('offsetWidth'),
                height: element.prop('offsetHeight'),
                top: boundingClientRect.top + ($window.pageYOffset || $document[0].body.scrollTop),
                left: boundingClientRect.left + ($window.pageXOffset || $document[0].body.scrollLeft)
            };
        }
    };
}]);








/**
 * @example
 * var x = $keypress.on('keypress', {'enter' : 'foo()'}, $scope);
 * ...
 * $keypress.off(x);
 */
Application.Interface.service('$keypress', ['keypressHelper', '$document', function(keypressHelper, $document){
    return {
        on : function(mode, actions, scope) {
            keypressHelper(mode, scope, $document, actions);
        },
        //todo -- see document.unbind() for format, use passed in scope
        //note - format: [elm, mode, handler, combinations];
        off:function(handle) {
            handle[0].unbind(handle[1], handle[2]);
        }
    }
}]);

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
        var handler = function (event) {
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
        };

        elm.bind(mode, handler);

        return [elm, mode, handler, combinations];
    };
}]);








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
                this.$scope.$destroy();

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

    //todo - more intelligent handling??
    Clotho.listen('serverAlert', function() {
        $scope.close('Another alert appeared');
        Clotho.say($scope.message);
    }, $scope);

}]);