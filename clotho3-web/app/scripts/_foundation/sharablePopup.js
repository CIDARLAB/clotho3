angular.module('clotho.clothoDirectives')
/**
 * @name sharable-popup
 *
 * @description Displays a popup showing the basic fields of an instance. Appended to Body.
 *
 * @note - need to create isolate scope so properties don't overlap... todo - try assigning to locals instead
 *
 * @attrs
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
			'placement="{{popup_placement}}" '+
			'>'+
			'</div>';

		return {
			restrict : 'EA',
			scope : {
				popupOpen : '=?' + prefix + 'Open'
			},
			controller : function ($scope, $element, $attrs) {

			},
			compile:function (tElement, tAttrs){
				var popupLinker = $compile(template);

				return function (scope, element, attrs, nullCtrl) {

					//short circuit link if no ID
					//todo - allow passing in of model
					if (angular.isUndefined(attrs[prefix + 'Id'])) {
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

						// Now set the calculated positioning.
						popup.css(popupPosition);
					}

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

					scope.sharable_id = "loading...";
					attrs.$observe(prefix + 'Id', function (val, oldval) {
						if (!!val && (!oldval || val != oldval)) {
							scope.sharable_id = val;
						}
					});

					attrs.$observe(prefix + 'Placement', function (val) {
						scope.popup_placement = angular.isDefined(val) ? val : '';
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

					scope.$watch('popupOpen', function (newval, oldval) {
						//hack - need to do this out of the $scope lifecycle because may need to trigger $digest
						setTimeout(function () {
							toggle(newval);
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
						removePopup();
						popup = null;
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
				sharableId: '=',
				placement: '@'
			},
			templateUrl: 'views/_foundation/sharableBasicFieldsPopup.html',
			link : function (scope, element, attrs, nullCtrl) {

				scope.$watch('sharableId', function ( val, oldval ) {
					if (!!val) {
						Clotho.get(val).then(function (retrievedSharable) {
							scope.fullSharable = retrievedSharable;
							scope.sharable = ClothoSchemas.pruneToBasicFields(retrievedSharable);
							scope.type = ClothoSchemas.determineInstanceType(retrievedSharable);
							scope.labelClass = 'label-' + ClothoSchemas.typeToColorClass(scope.type);

							if (ClothoSchemas.isSchema(retrievedSharable)) {
								scope.isSchema = true;
								ClothoSchemas.downloadSchemaDependencies(retrievedSharable).then(function (finalSchema) {
									scope.schema = finalSchema;
								});
							}
						});
					}
				});

				scope.$on('$destroy', function () {
				})
			}
		};
	});