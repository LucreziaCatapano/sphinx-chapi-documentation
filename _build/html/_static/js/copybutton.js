document.addEventListener('DOMContentLoaded', function() {
    const codeBlocks = document.querySelectorAll('.literal-block pre'); // Adjust selector if needed

    codeBlocks.forEach(function(codeBlock) {
        const copyButton = document.createElement('button');
        copyButton.classList.add('copy-button');
        copyButton.innerHTML = '<i class="fa fa-clipboard"></i>'; // Using Font Awesome icon
        copyButton.setAttribute('title', 'Copy to clipboard');

        const wrapper = document.createElement('div');
        wrapper.classList.add('code-block-wrapper'); // For positioning
        codeBlock.parentNode.insertBefore(wrapper, codeBlock);
        wrapper.appendChild(codeBlock);
        wrapper.appendChild(copyButton);

        const clipboard = new ClipboardJS(copyButton, {
            text: function(trigger) {
                return trigger.previousElementSibling.textContent; // Get text from <pre>
            }
        });

        clipboard.on('success', function(e) {
            // Optional: Provide feedback to the user
            const originalLabel = e.trigger.innerHTML;
            e.trigger.innerHTML = '<i class="fa fa-check" style="color: green;"></i> Copied!';
            setTimeout(function() {
                e.trigger.innerHTML = originalLabel;
            }, 2000);
            e.clearSelection();
        });

        clipboard.on('error', function(e) {
            // Optional: Handle errors
            e.trigger.innerHTML = 'Error';
        });
    });
});