//note - polyfill until migrate to $modal
//ignoring breaking up controllers etc. since will get rid of this

angular.module('clotho.interface').provider("$dialog", function(){

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
						templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
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

				login : function() {
					return new Dialog({
						backdrop: true,
						backdropFade: true,
						keyboard: true,
						backdropClick: true,
						templateUrl:  'views/_interface/ui-custom/dialogLogin.html',
						controller: 'DialogLoginController'
					});
				},

				serverAlert: function(message) {
					return new Dialog({
						backdrop: true,
						backdropFade: true,
						keyboard: true,
						backdropClick: true,
						templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
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
				},

				share : function(url) {
					return new Dialog({
						backdrop: true,
						backdropFade: true,
						keyboard: true,
						backdropClick: true,
						templateUrl:  'views/_interface/ui-custom/dialogShare.html',
						controller: 'DialogShareController',
						resolve: {
							model: function() {
								return {
									url : url
								}
							}
						}
					});
				},

				video : function(videoId, videoParams) {
					//todo: preserve aspect ratio
					angular.extend(videoParams, {width: "560"});

					return new Dialog({
						backdrop: true,
						backdropFade: true,
						keyboard: true,
						backdropClick: true,
						template:  '<div youtube="{{ videoId }}" params="videoParams""></div>',
						controller: 'VideoDialogController',
						resolve: {
							model: function() {
								return {
									videoId : videoId,
									videoParams : videoParams
								}
							}
						}
					})
				}
			};
		}];
});

angular.module('clotho.interface').controller('MessageBoxController', ['$scope', 'dialog', 'model', function($scope, dialog, model){
	$scope.title = model.title;
	$scope.message = model.message;
	$scope.buttons = model.buttons;
	$scope.close = function(res){
		dialog.close(res);
	};
}]);

angular.module('clotho.interface').controller('DialogLoginController', ['$scope', 'dialog', 'Clotho', function($scope, dialog, Clotho){
	$scope.close = function(res){
		dialog.close(res);
	};

	$scope.notification = {};
	$scope.cred = {username : "", password: ""};

	$scope.login = function() {
		Clotho.login($scope.cred.username, $scope.cred.password).then(function (result) {
			console.log('run login');
			if (!!result) {
				$scope.notification = {class : "alert-success", message: "Log in Success"};
				dialog.close($scope.cred.username);
			} else {
				$scope.notification = {class : "alert-danger", message: "Log in Error"};
				$scope.cred = {username : "", password: ""};
			}
		});
	};

}]);

angular.module('clotho.interface').controller('ServerAlertController', ['$scope', 'dialog', 'model', 'Clotho', function($scope, dialog, model, Clotho){
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

angular.module('clotho.interface').controller('DialogShareController', ['$scope', 'dialog', 'model', '$location', function($scope, dialog, model, $location){
	$scope.close = function(result){
		dialog.close(result);
	};

	$scope.customUrl = (model.url && model.url != '') ? model.url : false;

	$scope.social = [
		{
			"name" : "facebook",
			"prefix" : "http://www.facebook.com/sharer.php?u="
		},
		{
			"name" : "google",
			"prefix" : "https://plus.google.com/share?url="
		},
		{
			"name" : "twitter",
			"prefix" : "http://twitter.com/share?url="
		},
		{
			"name" : "linkedin",
			"prefix" : "http://www.linkedin.com/shareArticle?mini=true&url="
		},
		{
			"name" : "digg",
			"prefix" : "http://www.digg.com/submit?url="
		},
		{
			"name" : "reddit",
			"prefix" : "http://reddit.com/submit?url="
		},
		{
			"name" : "email",
			"prefix" : "mailto:?Body="
		}
	];

	$scope.share = function (site) {
		var url = $scope.customUrl ? $scope.customUrl : site.prefix + $location.absUrl();

		$scope.close();

		window.open(url, (site.name == 'email' ? '_self' : "_blank") );
	}

}]);


angular.module('clotho.interface').controller('VideoDialogController', ['$scope', 'dialog', 'model', function($scope, dialog, model){
	$scope.videoId = model.videoId;
	$scope.videoParams = model.videoParams;
	$scope.close = function(res){
		dialog.close(res);
	};
}]);