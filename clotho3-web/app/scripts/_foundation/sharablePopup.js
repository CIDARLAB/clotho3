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
 * sharablePopupTrigger
 * sharablePopupStartOpen //todo
 *
 * @usage
 *
 *
 *
 */

	//todo - incorporate with tokens for popover style

	.directive('sharablePopup', function ($animate, $window, $document, $compile, $timeout, Clotho, ClothoSchemas, hotkeys) {

		var bodyElement = $document.find( 'body' );

		var bodyOffset = function (element) {
			var boundingClientRect = element[0].getBoundingClientRect();
			return {
				width: boundingClientRect.width || element.prop('offsetWidth'),
				height: boundingClientRect.height || element.prop('offsetHeight'),
				top: boundingClientRect.top + ($window.pageYOffset || $document[0].body.scrollTop || $document[0].documentElement.scrollTop),
				left: boundingClientRect.left + ($window.pageXOffset || $document[0].body.scrollLeft  || $document[0].documentElement.scrollLeft)
			};
		};

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
			'sharable="tt_sharable" '+
			'placement="tt_placement" '+
			'>'+
			'</div>';

		return {
			restrict : 'EA',
			replace : true,
			scope : true,
			controller : function ($scope, $element, $attrs) {

			},
			link: function (scope, element, attrs, nullCtrl) {

				if (!attrs[prefix + 'Id']) {
					return;
				}

				var tooltip,
					tooltipLinker = $compile( template),
					triggers,
					hasRegisteredTriggers = false;

				function positionTooltip () {
					var position,
						ttWidth,
						ttHeight,
						ttPosition;
					// Get the position of the directive element.
					position = bodyOffset( element );

					// Get the height and width of the tooltip so we can center it.
					ttWidth = tooltip.prop( 'offsetWidth' );
					ttHeight = tooltip.prop( 'offsetHeight' );

					// Calculate the tooltip's top and left coordinates to center it with
					// this directive.
					switch ( scope.tt_placement ) {
						case 'right':
							ttPosition = {
								top: position.top + position.height / 2 - ttHeight / 2,
								left: position.left + position.width
							};
							break;
						case 'top':
							ttPosition = {
								top: position.top - ttHeight,
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
								top: position.top + position.height,
								left: position.left + position.width / 2 - ttWidth / 2
							};
							break;
					}

					ttPosition.top += 'px';
					ttPosition.left += 'px';

					// Now set the calculated positioning.
					tooltip.css( ttPosition );
				}

				function toggleTooltipBind () {
					if ( !scope.tt_isOpen ) {
						showTooltipBind();
					} else {
						hideTooltipBind();
					}
				}

				// Show the tooltip
				// double function to position correctly, see show return
				function showTooltipBind() {
					show()();
				}

				function hideTooltipBind () {
					scope.$apply(function () {
						hide();
					});
				}

				// Show the tooltip popup element.
				function show() {

					createTooltip();

					// Set the initial positioning.
					tooltip.css({ top: 0, left: 0, display: 'block' });

					//not visible anyway
					$animate.enter(tooltip, bodyElement , angular.element(bodyElement[0].lastChild), function () {});

					positionTooltip();

					// And show the tooltip.
					scope.tt_isOpen = true;
					scope.$digest(); // digest required as $apply is not called

					hotkeys.add('esc', hideTooltipBind);

					// Return positioning function as promise callback for correct
					// positioning after draw.
					return positionTooltip;
				}

				// Hide the tooltip popup element.
				function hide() {
					// First things first: we don't show it anymore.
					scope.tt_isOpen = false;
					removeTooltip();
					hotkeys.del('esc');
				}

				function createTooltip() {
					// There can only be one tooltip element per directive shown at once.
					if (tooltip) {
						removeTooltip();
					}
					tooltip = tooltipLinker(scope);

					// Get contents rendered into the tooltip
					scope.$digest();
				}

				function removeTooltip() {
					if (tooltip) {
						$animate.leave(tooltip, function () {
							tooltip = null;
						});
					}
				}

				/* WATCHERS */

				element.css({cursor : 'pointer'});

				scope.tt_sharable = "loading...";
				attrs.$observe( prefix+'Id', function ( val, oldval ) {
					//verify
					if (!!val && (!oldval || val != oldval)) {
						console.log(val);
						Clotho.get(val).then(function (retrievedSharable) {
							scope.tt_sharable = ClothoSchemas.pruneToBasicFields(retrievedSharable);
						});
					}
				});

				attrs.$observe( prefix+'Title', function ( val ) {
					scope.tt_title = val;
				});

				attrs.$observe( prefix+'Placement', function ( val ) {
					scope.tt_placement = angular.isDefined( val ) ? val : 'bottom';
				});

				var unregisterTriggers = function() {
					if (hasRegisteredTriggers) {
						element.unbind( triggers.show, showTooltipBind );
						element.unbind( triggers.hide, hideTooltipBind );
					}
				};

				attrs.$observe( prefix+'Trigger', function ( val ) {
					unregisterTriggers();

					triggers = getTriggers( val );

					if ( triggers.show === triggers.hide ) {
						element.bind( triggers.show, toggleTooltipBind );
					} else {
						element.bind( triggers.show, showTooltipBind );
						element.bind( triggers.hide, hideTooltipBind );
					}

					hasRegisteredTriggers = true;
				});

				//start closed unless specify attribute
				scope.tt_isOpen = angular.isDefined( attrs[prefix + 'StartOpen'] );
				if (scope.tt_isOpen) {
					//todo
				}

				// if a tooltip is attached to <body> we need to remove it on
				// location change as its parent scope will probably not be destroyed
				// by the change.
				scope.$on('$locationChangeSuccess', function closeTooltipOnLocationChangeSuccess () {
					if ( scope.tt_isOpen ) {
						hide();
					}
				});

				// Make sure tooltip is destroyed and removed.
				scope.$on('$destroy', function onDestroyTooltip() {
					unregisterTriggers();
					removeTooltip();
				});
			}
		}
	})
	.directive('sharablePopupInner', function () {
		return {
			restrict: 'EA',
			replace: true,
			scope: {
				sharable: '=',
				placement: '='
			},
			templateUrl: 'views/_foundation/sharableBasicFieldsPopup.html'
		};
	});