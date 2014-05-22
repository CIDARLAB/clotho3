angular.module('clotho.clothoDirectives')
/**
 * @name sharable-popup
 *
 * @description Displays a popup showing the basic fields of an instance. Appended to Body. Either pass Model or ID of sharable, model gets priority (should only use pruned fields or may be very large)
 *
 * @note - need to create isolate scope so properties don't overlap
 *
 * @attrs
 * sharablePopupModel
 * sharablePopupId
 * sharablePopupPlacement
 * sharablePopupTrigger (none | click | mouseenter | focus)
 * sharablePopupOpen
 *
 * @usage
 *
 *
 *
 */

	.directive('sharablePopup', function ($animate, $window, $document, $compile, $timeout, $parse, Clotho, ClothoSchemas, hotkeys) {

		var bodyElement = $document.find( 'body' );

		var bodyOffset = function (element) {
			var boundingClientRect = element[0].getBoundingClientRect();
			return {
				width: boundingClientRect.width || element.prop('offsetWidth'),
				height: boundingClientRect.height || element.prop('offsetHeight'),
				top: boundingClientRect.top + ($window.pageYOffset || $document[0].documentElement.scrollTop),
				left: boundingClientRect.left + ($window.pageXOffset || $document[0].documentElement.scrollLeft)
			};
		};

		// Default hide triggers for each show trigger
		var triggerMap = {
			'mouseenter': 'mouseleave',
			'click': 'click',
			'focus': 'blur'
		};

		function getTriggers ( trigger ) {
			var show = trigger || 'click';
			var hide = triggerMap[show] || show;
			return {
				show: show,
				hide: hide
			};
		}

		var prefix = 'sharablePopup';

		var template =
			'<div sharable-popup-inner '+
			'sharable-id="sharable_id" '+
			'sharable-model="passedModel" '+
			'placement="{{popup_placement}}" '+
			'reposition="repositionFunction()"' +
			'>'+
			'</div>';

		return {
			restrict : 'EA',
			scope : {
				passedModel : '=?' + prefix + 'Model'
			},
			compile:function (tElement, tAttrs){
				var popupLinker = $compile(template);

				return function (scope, element, attrs, nullCtrl) {

					//short circuit link if no ID or model
					if (angular.isUndefined(attrs[prefix + 'Id']) && angular.isUndefined(attrs[prefix + 'Model'])) {
						return;
					}

					var popup,
						popupScope,
						triggers,
						hasRegisteredTriggers = false;

					function positionPopup() {
						var position,
							popupWidth,
							popupHeight,
							popupPosition;
						// Get the position of the directive element.
						position = bodyOffset(element);

						// Get the height and width of the popup so we can center it.
						popupWidth = popup.prop('offsetWidth');
						popupHeight = popup.prop('offsetHeight');

						// Calculate the popup's top and left coordinates to center it with
						// this directive.
						switch (scope.popup_placement) {
							//e.g. in command bar
							case 'topRight':
								popupPosition = {
									top: position.top + position.height / 2 - 55,
									left: position.left + position.width
								};
								break;
							case 'right':
								popupPosition = {
									top: position.top + position.height / 2 - popupHeight / 2,
									left: position.left + position.width
								};
								break;
							case 'top':
								popupPosition = {
									top: position.top - popupHeight,
									left: position.left + position.width / 2 - popupWidth / 2
								};
								break;
							case 'left':
								popupPosition = {
									top: position.top + position.height / 2 - popupHeight / 2,
									left: position.left - popupWidth
								};
								break;
							case 'bottom':
								popupPosition = {
									top: position.top + position.height,
									left: position.left + position.width / 2 - popupWidth / 2
								};
								break;
							//display just under
							default :
							{
								popupPosition = {
									top: position.top + position.height,
									left: position.left
								};
								break;
							}
						}

						popupPosition.top += 'px';
						popupPosition.left += 'px';

						//check the topRight class, add right if current
						popup.toggleClass('right', scope.popup_placement == 'topRight');

						// Now set the calculated positioning.
						popup.css(popupPosition);
					}

					//expose on the scope so can be passed to the inner function to trigger when it's content change
					scope.repositionFunction = function () {
						$timeout(function () {
							positionPopup();
						});
					};

					function togglePopupBind() {
						if (!scope.popupOpen) {
							showPopupBind();
						} else {
							hidePopupBind();
						}
					}

					// Show the popup
					// double function to position correctly, see show return
					function showPopupBind() {
						show()();
					}

					function hidePopupBind() {
						scope.$apply(function () {
							hide();
						});
					}

					//NB toggleBind above... this is for internal use, can pass value to force
					function toggle(value) {
						if (angular.isDefined(value)) {
							!!value ? show()() : hide();
						} else {
							scope.popupOpen ? show()() : hide();
						}
					}

					// Show the popup element.
					function show() {

						createPopup();

						// Set the initial positioning.
						popup.css({ top: 0, left: 0, display: 'block' });

						//not visible anyway
						$animate.enter(popup, bodyElement, angular.element(bodyElement[0].lastChild), angular.noop);

						positionPopup();

						// And show the popup.
						scope.popupOpen = true;

						//bind hotkey to close
						hotkeys.add('esc', hide);

						// Return positioning function as promise callback for correct
						// positioning after draw.
						return positionPopup;
					}

					// Hide the popup element.
					function hide() {

						if (scope.popupOpen) {
							scope.popupOpen = false;
						}

						if (popup) {
							removePopup();
						}

						hotkeys.del('esc');
					}

					function createPopup() {
						if (popup) {
							removePopup();
						}

						// Make sure to use a new child scope every time as watchers leak into scope.
						// If linked DOM is removed, watchers from that DOM isn't removed.
						// Store it for manual destruction later
						popupScope = scope.$new();
						popup = popupLinker(popupScope, function () {});

						// Get contents rendered into the popup to position properly
						popupScope.$digest(); //because compiling ourselves
					}

					function removePopup() {
						if (popup) {
							popup.remove();
							popup = null;
						}
						if (popupScope) {
							popupScope.$destroy();
							popupScope = null;
						}
					}

					/* WATCHERS */

					element.css({cursor: 'pointer'});

					attrs.$observe(prefix + 'Id', function (val, oldval) {
						if (!!val && (!oldval || val != oldval)) {
							scope.sharable_id = val;
						}
					});

					attrs.$observe(prefix + 'Placement', function (val) {
						scope.popup_placement = angular.isDefined(val) ? val : 'topRight';
					});

					var unregisterTriggers = function () {
						if (hasRegisteredTriggers) {
							element.unbind(triggers.show, showPopupBind);
							element.unbind(triggers.hide, hidePopupBind);
						}
					};

					attrs.$observe(prefix + 'Trigger', function (val) {
						unregisterTriggers();

						if (val == 'none') {
							return;
						}

						triggers = getTriggers(val);

						if (triggers.show === triggers.hide) {
							element.bind(triggers.show, togglePopupBind);
						} else {
							element.bind(triggers.show, showPopupBind);
							element.bind(triggers.hide, hidePopupBind);
						}

						hasRegisteredTriggers = true;
					});

					scope.$watch(function () {
						return attrs[prefix + 'Open'];
					}, function (val) {
						scope.popupOpen = scope.$eval(val);
						setTimeout(function () {
							toggle(scope.popupOpen);
						});
					});

					// if a popup is attached to <body> we need to remove it on
					// location change as its parent scope will probably not be destroyed
					// by the change.
					scope.$on('$locationChangeSuccess', function closePopupOnLocationChangeSuccess() {
						if (scope.popupOpen) {
							hide();
						}
					});

					// Make sure popup is destroyed and removed.
					scope.$on('$destroy', function onDestroyPopup() {
						unregisterTriggers();
						hide();
					});
				}
			}
		}
	})
	.directive('sharablePopupInner', function (Clotho, ClothoSchemas) {
		return {
			restrict: 'EA',
			replace: true,
			scope: {
				sharableId: '=?',
				sharableModel : '=?',
				placement: '@',
				reposition : '&'
			},
			templateUrl: 'views/_foundation/sharableBasicFieldsPopup.html',
			link : function (scope, element, attrs, nullCtrl) {

				//given a model (should be pruned), sets related scope variables
				function setSharable (model) {
					scope.sharable = model;
					scope.type = ClothoSchemas.determineSharableType(model);
					scope.iconClass = ClothoSchemas.determineSharableIcon(scope.type);
					scope.labelClass = 'label-' + ClothoSchemas.typeToColorClass(scope.type);

					if (ClothoSchemas.isSchema(model)) {
						scope.isSchema = true;
						//must have a proper schema to download schema dependencies
						Clotho.get(model.id, {mute: true})
						.then(function (fullModel) {
							ClothoSchemas.downloadSchemaDependencies(fullModel)
							.then(function (finalSchema) {
								scope.schema = finalSchema;
							});
						})

					}
				}

				scope.$watch('sharableModel', function ( val, oldval ) {
					if (!!val) {
						setSharable(val);
						//no need to reposition if we pass in the model, since the initial $digest will fill it
					}
				});

				scope.$watch('sharableId', function ( val, oldval ) {
					if (!!val && angular.isEmpty(scope.sharableModel)) {
						Clotho.get(val, {mute : true}).then(function (retrievedSharable) {
							scope.fullSharable = retrievedSharable;
							setSharable( ClothoSchemas.pruneToBasicFields(retrievedSharable) );
							//if we're getting it remotely, the size probably changed, so let's reposition it.
							scope.reposition();
						});
					}
				});

				//remove focus from autocomplete input... hard to hide popup.
				scope.toggleSchema = function (evt) {
					evt.preventDefault();
					scope.showingSchema = !scope.showingSchema;
					scope.reposition();
				};

				scope.edit = Clotho.edit;

				scope.$on('$destroy', function () {
				})
			}
		};
	});