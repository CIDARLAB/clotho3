/**
 * @name highlight
 * @example

 <pre class="prettyprint">
 &lt;label&gt;&lt;input type=&quot;checkbox&quot; ng-model=&quot;caseSensitive&quot;&gt; Case Sensitive?&lt;/label&gt;
 &lt;input placeholder=&quot;Custom Class Name&quot; value=&quot;ui-match&quot; ng-model=&quot;className&quot;&gt;
 &lt;input placeholder=&quot;Enter some text to highlight&quot; value=&quot;you&quot; ng-model=&quot;highlightText&quot;&gt;
 &lt;p ng-bind-html-unsafe=&quot;&#x27;Hello there, how are you today? I\'m fine thank you.&#x27; | highlight:highlightText:className:caseSensitive&quot;&gt;&lt;/p&gt;

 &lt;style&gt;
 .ui-match { background: yellow; }
 &lt;/style&gt;
 </pre>

 * @param text {string} haystack to search through
 * @param search {string} needle to search for
 * @param {string} className Class name to use. Defaults to 'ui-match'
 * @param [caseSensitive] {boolean} optional boolean to use case-sensitive searching
 */
angular.module('clotho.interface').filter('highlight', function () {
	return function (text, search, className, caseSensitive) {
		if (search || angular.isNumber(search)) {
			className = angular.isDefined(className) ? className : 'ui-match';
			text = text.toString();
			search = search.toString();
			if (caseSensitive) {
				return text.split(search).join('<span class="'+className+'">' + search + '</span>');
			} else {
				return text.replace(new RegExp(search, 'gi'), '<span class="'+className+'">$&</span>');
			}
		} else {
			return text;
		}
	};
});